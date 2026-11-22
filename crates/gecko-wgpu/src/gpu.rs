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
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
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
    _pad0: u32,
    _pad1: u32,
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

    #[allow(clippy::too_many_arguments)]
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
        let keys_bytes = (idx.keys.len() * std::mem::size_of::<crate::index::KeyU64Pair>()) as u64;
        if keys_bytes > self.max_storage_binding {
            return Err(format!(
                "Index keys buffer ({:.1} MB) exceeds WGPU max_storage_buffer_binding_size ({:.1} MB).",
                keys_bytes as f64 / (1024.0 * 1024.0),
                self.max_storage_binding as f64 / (1024.0 * 1024.0),
            ));
        }

        let buf_x = {
            let x_bytes = seq_x.seq.to_ascii_vec_padded_u32();
            self.device
                .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("seq_x"),
                    contents: &x_bytes,
                    usage: wgpu::BufferUsages::STORAGE,
                })
        };
        let buf_y = {
            let y_bytes = seq_y.seq.to_ascii_vec_padded_u32();
            self.device
                .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("seq_y"),
                    contents: &y_bytes,
                    usage: wgpu::BufferUsages::STORAGE,
                })
        };
        let buf_y_rc = {
            let y_rc_bytes = match seq_y_rc {
                Some(rc) => rc.seq.to_ascii_vec_padded_u32(),
                None => vec![0u8; 16],
            };
            self.device
                .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("seq_y_rc"),
                    contents: &y_rc_bytes,
                    usage: wgpu::BufferUsages::STORAGE,
                })
        };
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

        let query_bg = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("query_bg"),
            layout: &query_pipeline.get_bind_group_layout(0),
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: buf_q_params.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: buf_x.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: buf_y.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: buf_y_rc.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 4,
                    resource: buf_bucket_off.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 5,
                    resource: buf_keys.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 6,
                    resource: buf_vals.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 7,
                    resource: buf_hits.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 8,
                    resource: buf_hit_count.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 9,
                    resource: buf_dropped_bitmap.as_entire_binding(),
                },
            ],
        });

        let extend_bg = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("extend_bg"),
            layout: &extend_pipeline.get_bind_group_layout(0),
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: buf_e_params.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: buf_x.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: buf_y.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: buf_y_rc.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 4,
                    resource: buf_hits.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 5,
                    resource: buf_out_frags.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 6,
                    resource: buf_frag_count.as_entire_binding(),
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
        let len_y = seq_y.seq.len as u32;
        let flags = (u32::from(soft_mask))
            | (u32::from(compare_reverse) << 1)
            | (u32::from(idx.has_dropped) << 2);
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
                shift_bits: idx.shift_bits,
                tile_start,
                tile_end,
                max_hits,
                stride,
                flags,
                _pad0: 0,
                _pad1: 0,
                _pad2: 0,
            };
            self.queue
                .write_buffer(&buf_q_params, 0, bytemuck::bytes_of(&q_params));
            self.queue.write_buffer(&buf_hit_count, 0, &[0u8; 4]);

            let n_threads = (tile_end - tile_start).div_ceil(stride);
            let n_wg = n_threads.div_ceil(128).min(self.max_workgroups);

            let mut enc = self.device.create_command_encoder(&Default::default());
            {
                let mut pass = enc.begin_compute_pass(&Default::default());
                pass.set_pipeline(&query_pipeline);
                pass.set_bind_group(0, &query_bg, &[]);
                pass.dispatch_workgroups(n_wg, 1, 1);
            }
            self.queue.submit(Some(enc.finish()));

            let n_hits = self.read_u32_buffer(&buf_hit_count);
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
                    _pad0: 0,
                    _pad1: 0,
                };
                self.queue
                    .write_buffer(&buf_e_params, 0, bytemuck::bytes_of(&e_params));
                self.queue.write_buffer(&buf_frag_count, 0, &[0u8; 4]);

                let ext_wg = n_hits.div_ceil(128).min(self.max_workgroups);
                let mut enc = self.device.create_command_encoder(&Default::default());
                {
                    let mut pass = enc.begin_compute_pass(&Default::default());
                    pass.set_pipeline(&extend_pipeline);
                    pass.set_bind_group(0, &extend_bg, &[]);
                    pass.dispatch_workgroups(ext_wg, 1, 1);
                }
                self.queue.submit(Some(enc.finish()));

                let n_frags_raw = self.read_u32_buffer(&buf_frag_count);
                if n_frags_raw > max_frags && (tile_end - tile_start) > 1 {
                    let mid = tile_start + (tile_end - tile_start) / 2;
                    tile_stack.push((mid, tile_end));
                    tile_stack.push((tile_start, mid));
                    continue;
                }
                let n_frags = n_frags_raw.min(max_frags) as usize;
                if n_frags > 0 {
                    let tile_frags: Vec<FragOut> = self.read_vec_buffer(&buf_out_frags, n_frags);
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
}
