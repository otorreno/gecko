use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use flate2::Compression;
use flate2::write::ZlibEncoder;

use crate::error::{GeckoError, Result};
use crate::formats::frags::Fragment;
use crate::seq::Genome;

/// Color encoding mode for dotplot fragment segments.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ColorMode {
    /// Blue for forward strand (`f`), red for reverse/inverted strand (`r`).
    #[default]
    Strand,
    /// Continuous heatmap colored by HSP similarity percentage (`0..=100`).
    Similarity,
}

/// Matrix origin convention.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PlotOrigin {
    /// Matrix style: `(0, 0)` at top-left, main forward diagonal runs top-left -> bottom-right
    /// (matches GECKO-MGV).
    #[default]
    TopLeft,
    /// Cartesian dotplot style: `(0, 0)` at bottom-left, main forward diagonal runs
    /// bottom-left -> top-right.
    BottomLeft,
}

/// Configuration for rendering a 2D synteny matrix PNG.
#[derive(Debug, Clone)]
pub struct PlotOptions {
    /// Total PNG image width in pixels (including margins).
    pub width: u32,
    /// Total PNG image height in pixels (including margins).
    pub height: u32,
    /// Stroke thickness in pixels for each HSP line segment.
    pub line_width: u32,
    /// Color encoding scheme.
    pub color_mode: ColorMode,
    /// Coordinate origin (`TopLeft` matrix style or `BottomLeft` Cartesian style).
    pub origin: PlotOrigin,
    /// Render with a dark background theme.
    pub dark: bool,
    /// Optional chart title.
    pub title: Option<String>,
    /// Optional X-axis label.
    pub x_label: Option<String>,
    /// Optional Y-axis label.
    pub y_label: Option<String>,
    /// Optional X contig boundaries (cumulative lengths excluding the final total).
    pub x_boundaries: Vec<u64>,
    /// Optional Y contig boundaries (cumulative lengths excluding the final total).
    pub y_boundaries: Vec<u64>,
}

impl Default for PlotOptions {
    fn default() -> Self {
        Self {
            width: 1200,
            height: 1200,
            line_width: 2,
            color_mode: ColorMode::Strand,
            origin: PlotOrigin::TopLeft,
            dark: false,
            title: None,
            x_label: None,
            y_label: None,
            x_boundaries: Vec::new(),
            y_boundaries: Vec::new(),
        }
    }
}

impl PlotOptions {
    /// Populate axis labels and multi-FASTA contig boundary guides from loaded [`Genome`]s.
    pub fn with_genomes(mut self, query: &Genome, reference: &Genome) -> Self {
        if self.x_label.is_none() {
            self.x_label = Some(format!("{} ({})", query.name, format_bp(query.total_bases)));
        }
        if self.y_label.is_none() {
            self.y_label = Some(format!(
                "{} ({})",
                reference.name,
                format_bp(reference.total_bases)
            ));
        }
        if self.title.is_none() {
            self.title = Some(format!("{} vs {}", query.name, reference.name));
        }
        self.x_boundaries = query
            .records
            .iter()
            .skip(1)
            .map(|r| r.start_no_sep)
            .collect();
        self.y_boundaries = reference
            .records
            .iter()
            .skip(1)
            .map(|r| r.start_no_sep)
            .collect();
        self
    }
}

/// Render fragments into a 2D synteny matrix PNG file on disk.
pub fn write_dotplot_png_path(
    path: impl AsRef<Path>,
    x_total: u64,
    y_total: u64,
    fragments: &[Fragment],
    opts: &PlotOptions,
) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::with_capacity(1 << 20, file);
    write_dotplot_png(&mut w, x_total, y_total, fragments, opts)?;
    w.flush()?;
    Ok(())
}

/// Render fragments into a 2D synteny matrix PNG stream.
pub fn write_dotplot_png<W: Write>(
    w: &mut W,
    x_total: u64,
    y_total: u64,
    fragments: &[Fragment],
    opts: &PlotOptions,
) -> Result<()> {
    if opts.width < 200 || opts.height < 200 {
        return Err(GeckoError::InvalidParameter(
            "plot width and height must be at least 200 pixels".into(),
        ));
    }

    let mut img = Canvas::new(opts.width, opts.height, palette(opts.dark).bg);
    let pal = palette(opts.dark);

    let margin_left: u32 = 96;
    let margin_right: u32 = 36;
    let margin_top: u32 = 64;
    let margin_bottom: u32 = 72;

    let plot_w = opts
        .width
        .saturating_sub(margin_left + margin_right)
        .max(32);
    let plot_h = opts
        .height
        .saturating_sub(margin_top + margin_bottom)
        .max(32);
    let px0 = margin_left as i32;
    let py0 = margin_top as i32;
    let px1 = px0 + plot_w as i32 - 1;
    let py1 = py0 + plot_h as i32 - 1;

    // Fill plot panel background
    img.fill_rect(px0, py0, plot_w as i32, plot_h as i32, pal.panel_bg);

    let x_max = x_total.max(1) as f64;
    let y_max = y_total.max(1) as f64;

    let map_x = |x: u64| -> i32 {
        let frac = (x as f64 / x_max).clamp(0.0, 1.0);
        px0 + (frac * (plot_w - 1) as f64).round() as i32
    };
    let map_y = |y: u64| -> i32 {
        let frac = (y as f64 / y_max).clamp(0.0, 1.0);
        match opts.origin {
            PlotOrigin::TopLeft => py0 + (frac * (plot_h - 1) as f64).round() as i32,
            PlotOrigin::BottomLeft => py1 - (frac * (plot_h - 1) as f64).round() as i32,
        }
    };

    // Coordinate gridlines & axis tick labels (5 divisions)
    let ticks = 5usize;
    for i in 0..=ticks {
        let frac = i as f64 / ticks as f64;
        let gx = px0 + (frac * (plot_w - 1) as f64).round() as i32;
        let gy = match opts.origin {
            PlotOrigin::TopLeft => py0 + (frac * (plot_h - 1) as f64).round() as i32,
            PlotOrigin::BottomLeft => py1 - (frac * (plot_h - 1) as f64).round() as i32,
        };

        if i > 0 && i < ticks {
            img.draw_vline(gx, py0, py1, pal.grid);
            img.draw_hline(px0, px1, gy, pal.grid);
        }

        // X tick + label
        img.draw_vline(gx, py1 + 1, py1 + 5, pal.axis);
        let x_bp = (frac * x_max).round() as u64;
        let x_txt = format_bp(x_bp);
        let tw = text_width(&x_txt, 1);
        let tx = (gx - tw / 2).clamp(4, opts.width as i32 - tw - 4);
        img.draw_text(tx, py1 + 10, &x_txt, 1, pal.text_muted);

        // Y tick + label
        img.draw_hline(px0 - 5, px0 - 1, gy, pal.axis);
        let y_bp = (frac * y_max).round() as u64;
        let y_txt = format_bp(y_bp);
        let yw = text_width(&y_txt, 1);
        let ty = (gy - 3).clamp(4, opts.height as i32 - 10);
        img.draw_text((px0 - 9 - yw).max(4), ty, &y_txt, 1, pal.text_muted);
    }

    // Multi-FASTA contig boundary dashed dividers
    for &xb in &opts.x_boundaries {
        let bx = map_x(xb);
        img.draw_dashed_vline(bx, py0, py1, pal.boundary);
    }
    for &yb in &opts.y_boundaries {
        let by = map_y(yb);
        img.draw_dashed_hline(px0, px1, by, pal.boundary);
    }

    // Sort indices by fragment length ascending so long syntenic blocks paint cleanly on top.
    let mut order: Vec<usize> = (0..fragments.len()).collect();
    order.sort_unstable_by_key(|&i| fragments[i].length);

    let radius = (opts.line_width.saturating_sub(1) / 2) as i32;
    let thickness = opts.line_width.max(1) as i32;

    for idx in order {
        let f = &fragments[idx];
        let x0 = map_x(f.x_start);
        let x1 = map_x(f.x_end);
        let (y_start_fwd, y_end_fwd) = if f.strand == 'r' {
            (
                y_total.saturating_sub(f.y_start + 1),
                y_total.saturating_sub(f.y_end + 1),
            )
        } else {
            (f.y_start, f.y_end)
        };
        let y0 = map_y(y_start_fwd);
        let y1 = map_y(y_end_fwd);

        let color = match opts.color_mode {
            ColorMode::Strand => {
                if f.strand == 'f' {
                    pal.fwd
                } else {
                    pal.rev
                }
            }
            ColorMode::Similarity => similarity_color(f.similarity, opts.dark),
        };

        img.draw_thick_line_clipped(x0, y0, x1, y1, radius, thickness, color, px0, py0, px1, py1);
    }

    // Plot border box
    img.draw_hline(px0, px1, py0, pal.axis);
    img.draw_hline(px0, px1, py1, pal.axis);
    img.draw_vline(px0, py0, py1, pal.axis);
    img.draw_vline(px1, py0, py1, pal.axis);

    // Title & legend in top margin
    let title = opts
        .title
        .clone()
        .unwrap_or_else(|| "GECKO 2D Synteny Matrix".to_string());
    img.draw_text(px0, 14, &title, 2, pal.text);

    let fwd_n = fragments.iter().filter(|f| f.strand == 'f').count();
    let rev_n = fragments.len() - fwd_n;
    let subtitle = format!(
        "HSPs: {} (fwd {}, rev {}) | X: {} | Y: {}",
        fragments.len(),
        fwd_n,
        rev_n,
        format_bp(x_total),
        format_bp(y_total)
    );
    img.draw_text(px0, 38, &subtitle, 1, pal.text_muted);

    // Legend swatches on top-right
    match opts.color_mode {
        ColorMode::Strand => {
            let legend_y = 38;
            let rev_label = "REV (-)";
            let fwd_label = "FWD (+)";
            let rx = px1 - text_width(rev_label, 1);
            img.draw_text(rx, legend_y, rev_label, 1, pal.text_muted);
            img.fill_rect(rx - 14, legend_y, 10, 7, pal.rev);

            let fx = rx - 28 - text_width(fwd_label, 1);
            img.draw_text(fx, legend_y, fwd_label, 1, pal.text_muted);
            img.fill_rect(fx - 14, legend_y, 10, 7, pal.fwd);
        }
        ColorMode::Similarity => {
            let bar_w = 100i32;
            let bx = px1 - bar_w;
            let by = 38;
            for dx in 0..bar_w {
                let sim = (dx as f32 / (bar_w - 1) as f32) * 100.0;
                let c = similarity_color(sim, opts.dark);
                img.draw_vline(bx + dx, by, by + 7, c);
            }
            let lbl = "SIM 0..100%";
            img.draw_text(bx - text_width(lbl, 1) - 8, by, lbl, 1, pal.text_muted);
        }
    }

    // X-axis title
    let x_title = opts
        .x_label
        .clone()
        .unwrap_or_else(|| format!("Query X ({})", format_bp(x_total)));
    let xtw = text_width(&x_title, 1);
    img.draw_text(
        px0 + (plot_w as i32 - xtw) / 2,
        opts.height as i32 - 22,
        &x_title,
        1,
        pal.text,
    );

    // Y-axis title (rendered compactly at top-left above Y ticks so it stays legible)
    let y_title = opts
        .y_label
        .clone()
        .unwrap_or_else(|| format!("Ref Y ({})", format_bp(y_total)));
    img.draw_text(8, py0 - 14, &truncate_str(&y_title, 26), 1, pal.text);

    encode_rgb_png(w, img.width, img.height, &img.pixels)
}

#[derive(Clone, Copy)]
struct Palette {
    bg: [u8; 3],
    panel_bg: [u8; 3],
    grid: [u8; 3],
    boundary: [u8; 3],
    axis: [u8; 3],
    text: [u8; 3],
    text_muted: [u8; 3],
    fwd: [u8; 3],
    rev: [u8; 3],
}

fn palette(dark: bool) -> Palette {
    if dark {
        Palette {
            bg: [15, 23, 42],        // slate-900
            panel_bg: [20, 30, 52],  // slate-850
            grid: [38, 52, 76],      // subtle grid
            boundary: [71, 85, 105], // contig divider
            axis: [148, 163, 184],   // slate-400
            text: [241, 245, 249],   // slate-100
            text_muted: [148, 163, 184],
            fwd: [56, 189, 248],  // sky-400
            rev: [251, 113, 133], // rose-400
        }
    } else {
        Palette {
            bg: [255, 255, 255],
            panel_bg: [250, 251, 253],
            grid: [226, 232, 240], // slate-200
            boundary: [186, 198, 214],
            axis: [51, 65, 85], // slate-700
            text: [15, 23, 42], // slate-900
            text_muted: [71, 85, 105],
            fwd: [29, 78, 216], // royal blue
            rev: [220, 38, 38], // crimson red
        }
    }
}

fn similarity_color(sim: f32, dark: bool) -> [u8; 3] {
    let t = (sim / 100.0).clamp(0.0, 1.0);
    // 4-stop perceptual gradient: indigo -> teal -> amber -> crimson
    let stops: [[f32; 3]; 4] = if dark {
        [
            [99.0, 102.0, 241.0],
            [45.0, 212.0, 191.0],
            [251.0, 191.0, 36.0],
            [244.0, 63.0, 94.0],
        ]
    } else {
        [
            [67.0, 56.0, 202.0],
            [13.0, 148.0, 136.0],
            [217.0, 119.0, 6.0],
            [190.0, 18.0, 60.0],
        ]
    };
    let scaled = t * 3.0;
    let idx = (scaled.floor() as usize).min(2);
    let local = scaled - idx as f32;
    let a = stops[idx];
    let b = stops[idx + 1];
    [
        (a[0] + (b[0] - a[0]) * local).round() as u8,
        (a[1] + (b[1] - a[1]) * local).round() as u8,
        (a[2] + (b[2] - a[2]) * local).round() as u8,
    ]
}

struct Canvas {
    width: u32,
    height: u32,
    pixels: Vec<u8>,
}

impl Canvas {
    fn new(width: u32, height: u32, bg: [u8; 3]) -> Self {
        let mut pixels = vec![0u8; (width as usize) * (height as usize) * 3];
        for px in pixels.as_chunks_mut::<3>().0 {
            *px = bg;
        }
        Self {
            width,
            height,
            pixels,
        }
    }

    #[inline]
    fn set_pixel(&mut self, x: i32, y: i32, color: [u8; 3]) {
        if x >= 0 && y >= 0 && (x as u32) < self.width && (y as u32) < self.height {
            let off = ((y as usize) * (self.width as usize) + (x as usize)) * 3;
            self.pixels[off] = color[0];
            self.pixels[off + 1] = color[1];
            self.pixels[off + 2] = color[2];
        }
    }

    fn fill_rect(&mut self, x: i32, y: i32, w: i32, h: i32, color: [u8; 3]) {
        for dy in 0..h {
            for dx in 0..w {
                self.set_pixel(x + dx, y + dy, color);
            }
        }
    }

    fn draw_hline(&mut self, x0: i32, x1: i32, y: i32, color: [u8; 3]) {
        for x in x0.min(x1)..=x0.max(x1) {
            self.set_pixel(x, y, color);
        }
    }

    fn draw_vline(&mut self, x: i32, y0: i32, y1: i32, color: [u8; 3]) {
        for y in y0.min(y1)..=y0.max(y1) {
            self.set_pixel(x, y, color);
        }
    }

    fn draw_dashed_hline(&mut self, x0: i32, x1: i32, y: i32, color: [u8; 3]) {
        for x in x0.min(x1)..=x0.max(x1) {
            if ((x - x0) & 7) < 4 {
                self.set_pixel(x, y, color);
            }
        }
    }

    fn draw_dashed_vline(&mut self, x: i32, y0: i32, y1: i32, color: [u8; 3]) {
        for y in y0.min(y1)..=y0.max(y1) {
            if ((y - y0) & 7) < 4 {
                self.set_pixel(x, y, color);
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn draw_thick_line_clipped(
        &mut self,
        mut x0: i32,
        mut y0: i32,
        x1: i32,
        y1: i32,
        radius: i32,
        thickness: i32,
        color: [u8; 3],
        min_x: i32,
        min_y: i32,
        max_x: i32,
        max_y: i32,
    ) {
        let dx = (x1 - x0).abs();
        let dy = -(y1 - y0).abs();
        let sx = if x0 < x1 { 1 } else { -1 };
        let sy = if y0 < y1 { 1 } else { -1 };
        let mut err = dx + dy;

        loop {
            if thickness <= 1 {
                if x0 >= min_x && x0 <= max_x && y0 >= min_y && y0 <= max_y {
                    self.set_pixel(x0, y0, color);
                }
            } else {
                let extra = if thickness % 2 == 0 { 1 } else { 0 };
                for oy in -radius..=(radius + extra) {
                    let py = y0 + oy;
                    if py < min_y || py > max_y {
                        continue;
                    }
                    for ox in -radius..=(radius + extra) {
                        let px = x0 + ox;
                        if px >= min_x && px <= max_x {
                            self.set_pixel(px, py, color);
                        }
                    }
                }
            }
            if x0 == x1 && y0 == y1 {
                break;
            }
            let e2 = 2 * err;
            if e2 >= dy {
                err += dy;
                x0 += sx;
            }
            if e2 <= dx {
                err += dx;
                y0 += sy;
            }
        }
    }

    fn draw_text(&mut self, x: i32, y: i32, text: &str, scale: i32, color: [u8; 3]) {
        let mut cx = x;
        for ch in text.chars() {
            let glyph = glyph_5x7(ch);
            for (row, &bits) in glyph.iter().enumerate() {
                for col in 0..5 {
                    if (bits & (1 << (4 - col))) != 0 {
                        self.fill_rect(
                            cx + col * scale,
                            y + (row as i32) * scale,
                            scale,
                            scale,
                            color,
                        );
                    }
                }
            }
            cx += 6 * scale;
        }
    }
}

fn text_width(text: &str, scale: i32) -> i32 {
    if text.is_empty() {
        0
    } else {
        (text.chars().count() as i32 * 6 - 1) * scale
    }
}

fn truncate_str(s: &str, max_chars: usize) -> String {
    if s.chars().count() <= max_chars {
        s.to_string()
    } else {
        let mut out: String = s.chars().take(max_chars.saturating_sub(2)).collect();
        out.push_str("..");
        out
    }
}

fn format_bp(bp: u64) -> String {
    if bp >= 1_000_000_000 {
        format!("{:.2}Gb", bp as f64 / 1e9)
    } else if bp >= 1_000_000 {
        format!("{:.2}Mb", bp as f64 / 1e6)
    } else if bp >= 1_000 {
        format!("{:.1}kb", bp as f64 / 1e3)
    } else {
        format!("{bp}bp")
    }
}

/// Minimal dependency-free 8-bit RGB PNG encoder using `flate2` (already in workspace).
fn encode_rgb_png<W: Write>(w: &mut W, width: u32, height: u32, rgb: &[u8]) -> Result<()> {
    // PNG magic signature
    w.write_all(&[137, 80, 78, 71, 13, 10, 26, 10])?;

    // IHDR chunk
    let mut ihdr = [0u8; 13];
    ihdr[0..4].copy_from_slice(&width.to_be_bytes());
    ihdr[4..8].copy_from_slice(&height.to_be_bytes());
    ihdr[8] = 8; // bit depth
    ihdr[9] = 2; // color type: Truecolor RGB
    ihdr[10] = 0; // compression: deflate
    ihdr[11] = 0; // filter method
    ihdr[12] = 0; // interlace: none
    write_png_chunk(w, b"IHDR", &ihdr)?;

    // IDAT chunk: prepend filter byte 0 to each RGB scanline and zlib-compress
    let row_bytes = (width as usize) * 3;
    let mut enc = ZlibEncoder::new(Vec::new(), Compression::fast());
    for y in 0..(height as usize) {
        enc.write_all(&[0u8])?;
        let start = y * row_bytes;
        enc.write_all(&rgb[start..start + row_bytes])?;
    }
    let compressed = enc.finish()?;
    write_png_chunk(w, b"IDAT", &compressed)?;

    // IEND chunk
    write_png_chunk(w, b"IEND", &[])?;
    Ok(())
}

fn write_png_chunk<W: Write>(w: &mut W, chunk_type: &[u8; 4], data: &[u8]) -> Result<()> {
    w.write_all(&(data.len() as u32).to_be_bytes())?;
    w.write_all(chunk_type)?;
    w.write_all(data)?;
    let crc = png_crc32(chunk_type, data);
    w.write_all(&crc.to_be_bytes())?;
    Ok(())
}

fn png_crc32(chunk_type: &[u8; 4], data: &[u8]) -> u32 {
    let mut crc: u32 = 0xffff_ffff;
    for &b in chunk_type.iter().chain(data.iter()) {
        crc ^= b as u32;
        for _ in 0..8 {
            if (crc & 1) != 0 {
                crc = (crc >> 1) ^ 0xedb8_8320;
            } else {
                crc >>= 1;
            }
        }
    }
    !crc
}

/// Compact built-in 5x7 ASCII bitmap font for axis ticks, titles, and legend.
fn glyph_5x7(ch: char) -> [u8; 7] {
    match ch.to_ascii_uppercase() {
        '0' => [0x0e, 0x11, 0x13, 0x15, 0x19, 0x11, 0x0e],
        '1' => [0x04, 0x0c, 0x04, 0x04, 0x04, 0x04, 0x0e],
        '2' => [0x0e, 0x11, 0x01, 0x06, 0x08, 0x10, 0x1f],
        '3' => [0x1e, 0x01, 0x01, 0x0e, 0x01, 0x01, 0x1e],
        '4' => [0x02, 0x06, 0x0a, 0x12, 0x1f, 0x02, 0x02],
        '5' => [0x1f, 0x10, 0x1e, 0x01, 0x01, 0x11, 0x0e],
        '6' => [0x06, 0x08, 0x10, 0x1e, 0x11, 0x11, 0x0e],
        '7' => [0x1f, 0x01, 0x02, 0x04, 0x08, 0x08, 0x08],
        '8' => [0x0e, 0x11, 0x11, 0x0e, 0x11, 0x11, 0x0e],
        '9' => [0x0e, 0x11, 0x11, 0x0f, 0x01, 0x02, 0x0c],
        'A' => [0x0e, 0x11, 0x11, 0x1f, 0x11, 0x11, 0x11],
        'B' => [0x1e, 0x11, 0x11, 0x1e, 0x11, 0x11, 0x1e],
        'C' => [0x0e, 0x11, 0x10, 0x10, 0x10, 0x11, 0x0e],
        'D' => [0x1c, 0x12, 0x11, 0x11, 0x11, 0x12, 0x1c],
        'E' => [0x1f, 0x10, 0x10, 0x1e, 0x10, 0x10, 0x1f],
        'F' => [0x1f, 0x10, 0x10, 0x1e, 0x10, 0x10, 0x10],
        'G' => [0x0e, 0x11, 0x10, 0x17, 0x11, 0x11, 0x0f],
        'H' => [0x11, 0x11, 0x11, 0x1f, 0x11, 0x11, 0x11],
        'I' => [0x0e, 0x04, 0x04, 0x04, 0x04, 0x04, 0x0e],
        'J' => [0x07, 0x02, 0x02, 0x02, 0x12, 0x12, 0x0c],
        'K' => [0x11, 0x12, 0x14, 0x18, 0x14, 0x12, 0x11],
        'L' => [0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x1f],
        'M' => [0x11, 0x1b, 0x15, 0x15, 0x11, 0x11, 0x11],
        'N' => [0x11, 0x19, 0x15, 0x13, 0x11, 0x11, 0x11],
        'O' => [0x0e, 0x11, 0x11, 0x11, 0x11, 0x11, 0x0e],
        'P' => [0x1e, 0x11, 0x11, 0x1e, 0x10, 0x10, 0x10],
        'Q' => [0x0e, 0x11, 0x11, 0x11, 0x15, 0x12, 0x0d],
        'R' => [0x1e, 0x11, 0x11, 0x1e, 0x14, 0x12, 0x11],
        'S' => [0x0f, 0x10, 0x10, 0x0e, 0x01, 0x01, 0x1e],
        'T' => [0x1f, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04],
        'U' => [0x11, 0x11, 0x11, 0x11, 0x11, 0x11, 0x0e],
        'V' => [0x11, 0x11, 0x11, 0x11, 0x11, 0x0a, 0x04],
        'W' => [0x11, 0x11, 0x11, 0x15, 0x15, 0x1b, 0x11],
        'X' => [0x11, 0x11, 0x0a, 0x04, 0x0a, 0x11, 0x11],
        'Y' => [0x11, 0x11, 0x0a, 0x04, 0x04, 0x04, 0x04],
        'Z' => [0x1f, 0x01, 0x02, 0x04, 0x08, 0x10, 0x1f],
        '.' => [0x00, 0x00, 0x00, 0x00, 0x00, 0x0c, 0x0c],
        ',' => [0x00, 0x00, 0x00, 0x00, 0x04, 0x04, 0x08],
        ':' => [0x00, 0x0c, 0x0c, 0x00, 0x0c, 0x0c, 0x00],
        '-' => [0x00, 0x00, 0x00, 0x1f, 0x00, 0x00, 0x00],
        '+' => [0x00, 0x04, 0x04, 0x1f, 0x04, 0x04, 0x00],
        '_' => [0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x1f],
        '/' => [0x01, 0x02, 0x02, 0x04, 0x08, 0x08, 0x10],
        '(' => [0x02, 0x04, 0x08, 0x08, 0x08, 0x04, 0x02],
        ')' => [0x08, 0x04, 0x02, 0x02, 0x02, 0x04, 0x08],
        '|' => [0x04, 0x04, 0x04, 0x04, 0x04, 0x04, 0x04],
        '%' => [0x19, 0x19, 0x02, 0x04, 0x08, 0x13, 0x13],
        '=' => [0x00, 0x00, 0x1f, 0x00, 0x1f, 0x00, 0x00],
        '<' => [0x02, 0x04, 0x08, 0x10, 0x08, 0x04, 0x02],
        '>' => [0x08, 0x04, 0x02, 0x01, 0x02, 0x04, 0x08],
        '[' => [0x0e, 0x08, 0x08, 0x08, 0x08, 0x08, 0x0e],
        ']' => [0x0e, 0x02, 0x02, 0x02, 0x02, 0x02, 0x0e],
        '~' => [0x00, 0x00, 0x08, 0x15, 0x02, 0x00, 0x00],
        _ => [0x00; 7],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encodes_valid_png_with_ihdr_idat_iend() {
        let frags = vec![
            Fragment {
                diag: 0,
                x_start: 100,
                y_start: 100,
                x_end: 800,
                y_end: 800,
                length: 701,
                ident: 690,
                score: 2716,
                similarity: 96.8,
                seq_x: 0,
                seq_y: 0,
                block: 0,
                strand: 'f',
            },
            Fragment {
                diag: 0,
                x_start: 200,
                y_start: 100,
                x_end: 500,
                y_end: 400,
                length: 301,
                ident: 290,
                score: 1116,
                similarity: 92.7,
                seq_x: 0,
                seq_y: 0,
                block: 0,
                strand: 'r',
            },
        ];

        let mut png = Vec::new();
        let opts = PlotOptions {
            width: 400,
            height: 400,
            ..Default::default()
        };
        write_dotplot_png(&mut png, 1000, 1000, &frags, &opts).unwrap();

        assert_eq!(&png[0..8], &[137, 80, 78, 71, 13, 10, 26, 10]);
        assert_eq!(&png[12..16], b"IHDR");
        assert_eq!(u32::from_be_bytes(png[16..20].try_into().unwrap()), 400);
        assert_eq!(u32::from_be_bytes(png[20..24].try_into().unwrap()), 400);
        assert!(png.windows(4).any(|w| w == b"IDAT"));
        assert_eq!(&png[png.len() - 8..png.len() - 4], b"IEND");
    }
}
