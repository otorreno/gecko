//! WGPU adapter discovery, device setup, and compute pipeline dispatch for
//! `query.wgsl` and `extend.wgsl`.

use bytemuck::{Pod, Zeroable};
use gecko_core::extend::Scoring;
use gecko_core::seq::Genome;
use wgpu::util::DeviceExt;

use crate::index::WgpuRefIndex;

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct QueryParams {
    len_x: u32,
    len_y: u32,
    wl: u32,
    shift_bits: u32,
    tile_start: u32,
    tile_end: u32,
    max_hits: u32,
    stride: u32,
    flags: u32,
    inv_y_word_off: u32,
    _pad0: u32,
    _pad1: u32,
}

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct ExtendParams {
    len_x: u32,
    len_y: u32,
    wl: u32,
    min_len: u32,
    sim_th_x100: u32,
    n_hits: u32,
    max_frags: u32,
    match_score: i32,
    mismatch_penalty: i32,
    n_penalty: i32,
    sp_x_off: u32,
    sp_y_off: u32,
    sp_yrc_off: u32,
    sep_x_off: u32,
    sep_x_len: u32,
    sep_y_off: u32,
    sep_y_len: u32,
    sep_yrc_off: u32,
    sep_yrc_len: u32,
    iupac_x_off: u32,
    iupac_x_len: u32,
    iupac_y_off: u32,
    iupac_y_len: u32,
    iupac_yrc_off: u32,
    iupac_yrc_len: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
}

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
pub struct HitRec {
    pub pos_x: u32,
    pub pos_y: u32,
    pub strand: u32,
    pub _pad: u32,
}

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
pub struct FragOut {
    pub x_start: u32,
    pub y_start: u32,
    pub x_end: u32,
    pub y_end: u32,
    pub ident: u32,
    pub score: u32,
    pub strand: u32,
    pub _pad: u32,
}

pub struct GpuRunStats {
    pub leftmax_hits: u64,
    pub extended_hits: u64,
    pub raw_frags: Vec<FragOut>,
}

pub struct WgpuContext {
    pub adapter_name: String,
    pub adapter_backend: String,
    pub is_software: bool,
    device: wgpu::Device,
    queue: wgpu::Queue,
    max_storage_binding: u64,
    max_workgroups: u32,
}

impl WgpuContext {
    pub fn init(allow_software: bool) -> Result<Self, String> {
        let backends = wgpu::Backends::VULKAN | wgpu::Backends::METAL | wgpu::Backends::DX12;
        let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor {
            backends,
            ..Default::default()
        });

        let adapters = instance.enumerate_adapters(backends);

        let mut best: Option<(u8, wgpu::Adapter)> = None;
        for a in adapters {
            let info = a.get_info();
            let rank = match info.device_type {
                wgpu::DeviceType::DiscreteGpu => 3,
                wgpu::DeviceType::IntegratedGpu => 2,
                wgpu::DeviceType::VirtualGpu => 1,
                wgpu::DeviceType::Cpu | wgpu::DeviceType::Other => {
                    if allow_software {
                        0
                    } else {
                        continue;
                    }
                }
            };
            if best.as_ref().is_none_or(|(br, _)| rank > *br) {
                best = Some((rank, a));
            }
        }

        let adapter = best
            .map(|(_, a)| a)
            .ok_or_else(|| "No suitable WGPU GPU adapter found".to_string())?;
        let info = adapter.get_info();
        let is_software = matches!(info.device_type, wgpu::DeviceType::Cpu);

        let limits = adapter.limits();
        let max_storage_binding = limits.max_storage_buffer_binding_size as u64;
        let max_workgroups = limits.max_compute_workgroups_per_dimension;

        let (device, queue) = pollster::block_on(adapter.request_device(
            &wgpu::DeviceDescriptor {
                label: Some("gecko_wgpu_device"),
                required_features: wgpu::Features::empty(),
                required_limits: limits,
                memory_hints: wgpu::MemoryHints::Performance,
            },
            None,
        ))
        .map_err(|e| format!("Failed to create WGPU device: {e}"))?;

        Ok(Self {
            adapter_name: info.name,
            adapter_backend: format!("{:?}", info.backend),
            is_software,
            device,
            queue,
            max_storage_binding,
            max_workgroups,
        })
    }

    fn read_u32_buffer(&self, src: &wgpu::Buffer) -> u32 {
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("read_u32_staging"),
            size: 4,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let mut enc = self.device.create_command_encoder(&Default::default());
        enc.copy_buffer_to_buffer(src, 0, &staging, 0, 4);
        self.queue.submit(Some(enc.finish()));

        let slice = staging.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx.send(r);
        });
        self.device.poll(wgpu::Maintain::Wait);
        rx.recv().unwrap().unwrap();
        let data = slice.get_mapped_range();
        u32::from_le_bytes(data[0..4].try_into().unwrap())
    }

    fn read_vec_buffer<T: Pod>(&self, src: &wgpu::Buffer, count: usize) -> Vec<T> {
        if count == 0 {
            return Vec::new();
        }
        let byte_len = (count * std::mem::size_of::<T>()) as u64;
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("read_vec_staging"),
            size: byte_len,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let mut enc = self.device.create_command_encoder(&Default::default());
        enc.copy_buffer_to_buffer(src, 0, &staging, 0, byte_len);
        self.queue.submit(Some(enc.finish()));

        let slice = staging.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx.send(r);
        });
        self.device.poll(wgpu::Maintain::Wait);
        rx.recv().unwrap().unwrap();
        let data = slice.get_mapped_range();
        bytemuck::cast_slice(&data).to_vec()
    }
}

/// Persistent GPU-resident reference buffers, compiled WGSL compute pipelines, and reusable
/// scratch buffers for prepare-once, query-many execution (e.g. MetaGECKO read batches).
pub struct WgpuPreparedRef {
    len_y: u32,
    shift_bits: u32,
    has_dropped: bool,
    buf_y: wgpu::Buffer,
    buf_y_rc: wgpu::Buffer,
    buf_bucket_off: wgpu::Buffer,
    buf_keys: wgpu::Buffer,
    buf_vals: wgpu::Buffer,
    buf_dropped_bitmap: wgpu::Buffer,
    inv_y: Option<Vec<u32>>,
    sp_y: Option<Vec<u32>>,
    sp_yrc: Option<Vec<u32>>,
    sep_y: Vec<u32>,
    sep_yrc: Vec<u32>,
    iupac_y: Vec<u32>,
    iupac_yrc: Vec<u32>,
    query_pipeline: wgpu::ComputePipeline,
    extend_pipeline: wgpu::ComputePipeline,
    max_hits: u32,
    max_frags: u32,
    buf_hits: wgpu::Buffer,
    buf_hit_count: wgpu::Buffer,
    buf_out_frags: wgpu::Buffer,
    buf_frag_count: wgpu::Buffer,
    buf_q_params: wgpu::Buffer,
    buf_e_params: wgpu::Buffer,
}

impl WgpuContext {
    /// Upload reference sequences + Hash-CSR index and compile WGSL pipelines once on the GPU.
    pub fn prepare_ref(
        &self,
        seq_y: &Genome,
        seq_y_rc: Option<&Genome>,
        idx: &WgpuRefIndex,
        soft_mask: bool,
    ) -> Result<WgpuPreparedRef, String> {
        let keys_bytes = (idx.keys.len() * std::mem::size_of::<crate::index::KeyU64Pair>()) as u64;
        if keys_bytes > self.max_storage_binding {
            return Err(format!(
                "Index keys buffer ({:.1} MB) exceeds WGPU max_storage_buffer_binding_size ({:.1} MB).",
                keys_bytes as f64 / (1024.0 * 1024.0),
                self.max_storage_binding as f64 / (1024.0 * 1024.0),
            ));
        }

        let buf_y = self
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("seq_y_bits"),
                contents: bytemuck::cast_slice(&seq_y.seq.bits),
                usage: wgpu::BufferUsages::STORAGE,
            });
        let dummy_bits = [0u64; 2];
        let buf_y_rc = self
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("seq_y_rc_bits"),
                contents: bytemuck::cast_slice(match seq_y_rc {
                    Some(rc) => rc.seq.bits.as_slice(),
                    None => &dummy_bits,
                }),
                usage: wgpu::BufferUsages::STORAGE,
            });

        let buf_bucket_off = self
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("bucket_off"),
                contents: bytemuck::cast_slice(&idx.bucket_off),
                usage: wgpu::BufferUsages::STORAGE,
            });
        let buf_keys = self
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("idx_keys"),
                contents: bytemuck::cast_slice(&idx.keys),
                usage: wgpu::BufferUsages::STORAGE,
            });
        let buf_vals = self
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("idx_vals"),
                contents: bytemuck::cast_slice(&idx.vals),
                usage: wgpu::BufferUsages::STORAGE,
            });
        let buf_dropped_bitmap =
            self.device
                .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("dropped_bitmap"),
                    contents: bytemuck::cast_slice(&idx.dropped_bitmap),
                    usage: wgpu::BufferUsages::STORAGE,
                });

        let query_shader = self
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("query_wgsl"),
                source: wgpu::ShaderSource::Wgsl(include_str!("shaders/query.wgsl").into()),
            });
        let extend_shader = self
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("extend_wgsl"),
                source: wgpu::ShaderSource::Wgsl(include_str!("shaders/extend.wgsl").into()),
            });

        let query_pipeline =
            self.device
                .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                    label: Some("query_pipeline"),
                    layout: None,
                    module: &query_shader,
                    entry_point: Some("main"),
                    compilation_options: Default::default(),
                    cache: None,
                });
        let extend_pipeline =
            self.device
                .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                    label: Some("extend_pipeline"),
                    layout: None,
                    module: &extend_shader,
                    entry_point: Some("main"),
                    compilation_options: Default::default(),
                    cache: None,
                });

        let max_hits = ((self.max_storage_binding / 16).min(8_000_000)) as u32;
        let max_frags = ((self.max_storage_binding / 32).min(2_000_000)) as u32;

        let buf_hits = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("tile_hits"),
            size: (max_hits as u64) * std::mem::size_of::<HitRec>() as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let buf_hit_count = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("hit_count"),
            size: 4,
            usage: wgpu::BufferUsages::STORAGE
                | wgpu::BufferUsages::COPY_SRC
                | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let buf_out_frags = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("out_frags"),
            size: (max_frags as u64) * std::mem::size_of::<FragOut>() as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let buf_frag_count = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("frag_count"),
            size: 4,
            usage: wgpu::BufferUsages::STORAGE
                | wgpu::BufferUsages::COPY_SRC
                | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let buf_q_params = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("query_params"),
            size: std::mem::size_of::<QueryParams>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let buf_e_params = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("extend_params"),
            size: std::mem::size_of::<ExtendParams>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let inv_y = seq_y.seq.invalid_seed_bitmap(soft_mask);
        let sp_y = seq_y.seq.special_bitmap();
        let sp_yrc = seq_y_rc.and_then(|rc| rc.seq.special_bitmap());
        let sep_y: Vec<u32> = seq_y.seq.sep_pos.iter().map(|&p| p as u32).collect();
        let sep_yrc: Vec<u32> = match seq_y_rc {
            Some(rc) => rc.seq.sep_pos.iter().map(|&p| p as u32).collect(),
            None => Vec::new(),
        };
        let encode_iupac = |iupac: &[(u64, u8)]| -> Vec<u32> {
            let mut out = Vec::with_capacity(iupac.len() * 2);
            for &(p, b) in iupac {
                out.push(p as u32);
                out.push(b as u32);
            }
            out
        };
        let iupac_y = encode_iupac(&seq_y.seq.iupac);
        let iupac_yrc = match seq_y_rc {
            Some(rc) => encode_iupac(&rc.seq.iupac),
            None => Vec::new(),
        };

        Ok(WgpuPreparedRef {
            len_y: seq_y.seq.len as u32,
            shift_bits: idx.shift_bits,
            has_dropped: idx.has_dropped,
            buf_y,
            buf_y_rc,
            buf_bucket_off,
            buf_keys,
            buf_vals,
            buf_dropped_bitmap,
            inv_y,
            sp_y,
            sp_yrc,
            sep_y,
            sep_yrc,
            iupac_y,
            iupac_yrc,
            query_pipeline,
            extend_pipeline,
            max_hits,
            max_frags,
            buf_hits,
            buf_hit_count,
            buf_out_frags,
            buf_frag_count,
            buf_q_params,
            buf_e_params,
        })
    }

    /// Query a genome or metagenomic read batch against a persistent GPU-resident `WgpuPreparedRef`.
    #[allow(clippy::too_many_arguments)]
    pub fn query_prepared_ref(
        &self,
        seq_x: &Genome,
        prep: &WgpuPreparedRef,
        wl: u32,
        min_len: u32,
        min_similarity: f32,
        soft_mask: bool,
        compare_reverse: bool,
        scoring: Scoring,
    ) -> Result<GpuRunStats, String> {
        const NONE_OFF: u32 = u32::MAX;

        // Upload only the query batch (0.25 bytes/bp — 4x smaller than ASCII).
        let buf_x = self
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("seq_x_bits"),
                contents: bytemuck::cast_slice(&seq_x.seq.bits),
                usage: wgpu::BufferUsages::STORAGE,
            });

        // Pack invalid-seed bitmaps (X then cached Y) into a single storage buffer for query.wgsl.
        let inv_x = seq_x.seq.invalid_seed_bitmap(soft_mask);
        let has_inv_x = inv_x.is_some();
        let has_inv_y = prep.inv_y.is_some();
        let mut invalid_bm_data: Vec<u32> = Vec::new();
        if let Some(bm) = inv_x {
            invalid_bm_data.extend_from_slice(&bm);
        }
        let inv_y_word_off = invalid_bm_data.len() as u32;
        if let Some(ref bm) = prep.inv_y {
            invalid_bm_data.extend_from_slice(bm);
        }
        if invalid_bm_data.is_empty() {
            invalid_bm_data.push(0);
        }
        let buf_invalid_bm = self
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("invalid_bm"),
                contents: bytemuck::cast_slice(&invalid_bm_data),
                usage: wgpu::BufferUsages::STORAGE,
            });

        // Pack special bitmaps (X, cached Y, cached Y_rc) into a single storage buffer for extend.wgsl.
        let mut special_bm_data: Vec<u32> = Vec::new();
        let sp_x_off = if let Some(bm) = seq_x.seq.special_bitmap() {
            let off = special_bm_data.len() as u32;
            special_bm_data.extend_from_slice(&bm);
            off
        } else {
            NONE_OFF
        };
        let sp_y_off = if let Some(ref bm) = prep.sp_y {
            let off = special_bm_data.len() as u32;
            special_bm_data.extend_from_slice(bm);
            off
        } else {
            NONE_OFF
        };
        let sp_yrc_off = if let Some(ref bm) = prep.sp_yrc {
            let off = special_bm_data.len() as u32;
            special_bm_data.extend_from_slice(bm);
            off
        } else {
            NONE_OFF
        };
        if special_bm_data.is_empty() {
            special_bm_data.push(0);
        }
        let buf_special_bm = self
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("special_bm"),
                contents: bytemuck::cast_slice(&special_bm_data),
                usage: wgpu::BufferUsages::STORAGE,
            });

        // Pack sparse SEQ_SEP and rare IUPAC metadata (X, cached Y, cached Y_rc).
        let mut meta_data: Vec<u32> = Vec::new();
        let sep_x_off = meta_data.len() as u32;
        meta_data.extend(seq_x.seq.sep_pos.iter().map(|&p| p as u32));
        let sep_x_len = seq_x.seq.sep_pos.len() as u32;

        let sep_y_off = meta_data.len() as u32;
        meta_data.extend_from_slice(&prep.sep_y);
        let sep_y_len = prep.sep_y.len() as u32;

        let sep_yrc_off = meta_data.len() as u32;
        meta_data.extend_from_slice(&prep.sep_yrc);
        let sep_yrc_len = prep.sep_yrc.len() as u32;

        let iupac_x_off = meta_data.len() as u32;
        for &(p, b) in &seq_x.seq.iupac {
            meta_data.push(p as u32);
            meta_data.push(b as u32);
        }
        let iupac_x_len = seq_x.seq.iupac.len() as u32;

        let iupac_y_off = meta_data.len() as u32;
        meta_data.extend_from_slice(&prep.iupac_y);
        let iupac_y_len = (prep.iupac_y.len() / 2) as u32;

        let iupac_yrc_off = meta_data.len() as u32;
        meta_data.extend_from_slice(&prep.iupac_yrc);
        let iupac_yrc_len = (prep.iupac_yrc.len() / 2) as u32;

        if meta_data.is_empty() {
            meta_data.push(0);
        }
        let buf_meta = self
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("meta_buf"),
                contents: bytemuck::cast_slice(&meta_data),
                usage: wgpu::BufferUsages::STORAGE,
            });

        let query_bg = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("query_bg"),
            layout: &prep.query_pipeline.get_bind_group_layout(0),
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: prep.buf_q_params.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: buf_x.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: prep.buf_y.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: buf_invalid_bm.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 4,
                    resource: prep.buf_bucket_off.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 5,
                    resource: prep.buf_keys.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 6,
                    resource: prep.buf_vals.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 7,
                    resource: prep.buf_hits.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 8,
                    resource: prep.buf_hit_count.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 9,
                    resource: prep.buf_dropped_bitmap.as_entire_binding(),
                },
            ],
        });

        let extend_bg = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("extend_bg"),
            layout: &prep.extend_pipeline.get_bind_group_layout(0),
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: prep.buf_e_params.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: buf_x.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: prep.buf_y.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: prep.buf_y_rc.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 4,
                    resource: buf_special_bm.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 5,
                    resource: buf_meta.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 6,
                    resource: prep.buf_hits.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 7,
                    resource: prep.buf_out_frags.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 8,
                    resource: prep.buf_frag_count.as_entire_binding(),
                },
            ],
        });

        let mut init_tile_bp: u32 = if wl <= 14 {
            1024 * 1024
        } else if wl <= 19 {
            4 * 1024 * 1024
        } else {
            8 * 1024 * 1024
        };
        if !soft_mask && init_tile_bp > 1024 * 1024 {
            init_tile_bp >>= 1;
        }
        let stride: u32 = 16;
        let len_x = seq_x.seq.len as u32;
        let len_y = prep.len_y;
        let max_hits = prep.max_hits;
        let max_frags = prep.max_frags;
        let flags = u32::from(compare_reverse)
            | (u32::from(prep.has_dropped) << 1)
            | (u32::from(has_inv_x) << 2)
            | (u32::from(has_inv_y) << 3);
        let sim_th_x100 = (min_similarity * 100.0).round().max(0.0) as u32;

        let mut raw_frags = Vec::new();
        let mut tot_leftmax: u64 = 0;
        let mut tot_extended: u64 = 0;

        let mut tile_stack: Vec<(u32, u32)> = Vec::new();
        let mut ts: u32 = 0;
        while ts < len_x {
            let te = ts.saturating_add(init_tile_bp).min(len_x);
            tile_stack.push((ts, te));
            ts = te;
        }
        tile_stack.reverse();

        while let Some((tile_start, tile_end)) = tile_stack.pop() {
            let q_params = QueryParams {
                len_x,
                len_y,
                wl,
                shift_bits: prep.shift_bits,
                tile_start,
                tile_end,
                max_hits,
                stride,
                flags,
                inv_y_word_off,
                _pad0: 0,
                _pad1: 0,
            };
            self.queue
                .write_buffer(&prep.buf_q_params, 0, bytemuck::bytes_of(&q_params));
            self.queue.write_buffer(&prep.buf_hit_count, 0, &[0u8; 4]);

            let n_threads = (tile_end - tile_start).div_ceil(stride);
            let n_wg = n_threads.div_ceil(128).min(self.max_workgroups);

            let mut enc = self.device.create_command_encoder(&Default::default());
            {
                let mut pass = enc.begin_compute_pass(&Default::default());
                pass.set_pipeline(&prep.query_pipeline);
                pass.set_bind_group(0, &query_bg, &[]);
                pass.dispatch_workgroups(n_wg, 1, 1);
            }
            self.queue.submit(Some(enc.finish()));

            let n_hits = self.read_u32_buffer(&prep.buf_hit_count);
            if n_hits > max_hits && (tile_end - tile_start) > 1 {
                let mid = tile_start + (tile_end - tile_start) / 2;
                tile_stack.push((mid, tile_end));
                tile_stack.push((tile_start, mid));
                continue;
            }
            let n_hits = n_hits.min(max_hits);

            if n_hits > 0 {
                let e_params = ExtendParams {
                    len_x,
                    len_y,
                    wl,
                    min_len,
                    sim_th_x100,
                    n_hits,
                    max_frags,
                    match_score: scoring.match_score as i32,
                    mismatch_penalty: scoring.mismatch_penalty as i32,
                    n_penalty: scoring.n_penalty as i32,
                    sp_x_off,
                    sp_y_off,
                    sp_yrc_off,
                    sep_x_off,
                    sep_x_len,
                    sep_y_off,
                    sep_y_len,
                    sep_yrc_off,
                    sep_yrc_len,
                    iupac_x_off,
                    iupac_x_len,
                    iupac_y_off,
                    iupac_y_len,
                    iupac_yrc_off,
                    iupac_yrc_len,
                    _pad0: 0,
                    _pad1: 0,
                    _pad2: 0,
                };
                self.queue
                    .write_buffer(&prep.buf_e_params, 0, bytemuck::bytes_of(&e_params));
                self.queue.write_buffer(&prep.buf_frag_count, 0, &[0u8; 4]);

                let ext_wg = n_hits.div_ceil(128).min(self.max_workgroups);
                let mut enc = self.device.create_command_encoder(&Default::default());
                {
                    let mut pass = enc.begin_compute_pass(&Default::default());
                    pass.set_pipeline(&prep.extend_pipeline);
                    pass.set_bind_group(0, &extend_bg, &[]);
                    pass.dispatch_workgroups(ext_wg, 1, 1);
                }
                self.queue.submit(Some(enc.finish()));

                let n_frags_raw = self.read_u32_buffer(&prep.buf_frag_count);
                if n_frags_raw > max_frags && (tile_end - tile_start) > 1 {
                    let mid = tile_start + (tile_end - tile_start) / 2;
                    tile_stack.push((mid, tile_end));
                    tile_stack.push((tile_start, mid));
                    continue;
                }
                let n_frags = n_frags_raw.min(max_frags) as usize;
                if n_frags > 0 {
                    let tile_frags: Vec<FragOut> =
                        self.read_vec_buffer(&prep.buf_out_frags, n_frags);
                    raw_frags.extend_from_slice(&tile_frags);
                }
            }

            tot_leftmax += n_hits as u64;
            tot_extended += n_hits as u64;
        }

        Ok(GpuRunStats {
            leftmax_hits: tot_leftmax,
            extended_hits: tot_extended,
            raw_frags,
        })
    }

    #[allow(clippy::too_many_arguments, dead_code)]
    pub fn run_query_and_extend(
        &self,
        seq_x: &Genome,
        seq_y: &Genome,
        seq_y_rc: Option<&Genome>,
        idx: &WgpuRefIndex,
        wl: u32,
        min_len: u32,
        min_similarity: f32,
        soft_mask: bool,
        compare_reverse: bool,
        scoring: Scoring,
    ) -> Result<GpuRunStats, String> {
        let prep = self.prepare_ref(seq_y, seq_y_rc, idx, soft_mask)?;
        self.query_prepared_ref(
            seq_x,
            &prep,
            wl,
            min_len,
            min_similarity,
            soft_mask,
            compare_reverse,
            scoring,
        )
    }
}
