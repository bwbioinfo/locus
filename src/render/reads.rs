use ratatui::{
    buffer::Buffer,
    layout::Rect,
    style::{Color, Modifier, Style},
    widgets::Widget,
};

use crate::cache::{CigarOp, PileupRow, RenderRead, Strand};

use super::ViewTransform;

/// Show individual bases when the view is this many bp per column or narrower.
const BASE_RENDER_THRESHOLD: f64 = 5.0;

pub struct ReadsTrack<'a> {
    pub reads: &'a [RenderRead],
    pub rows: &'a [PileupRow],
    pub transform: ViewTransform,
    #[allow(dead_code)]
    pub show_names: bool,
}

impl<'a> Widget for ReadsTrack<'a> {
    fn render(self, area: Rect, buf: &mut Buffer) {
        if area.height == 0 || area.width == 0 {
            return;
        }

        let show_bases = self.transform.bp_per_col() <= BASE_RENDER_THRESHOLD;

        for (row_idx, row) in self.rows.iter().enumerate() {
            if row_idx >= area.height as usize {
                break;
            }
            let y = area.y + row_idx as u16;

            for &read_idx in row {
                let Some(read) = self.reads.get(read_idx) else {
                    continue;
                };
                if show_bases {
                    render_bases(read, y, area, &self.transform, buf);
                } else {
                    render_arrows(read, y, area, &self.transform, buf);
                }
            }
        }
    }
}

// ─── Base-level rendering (zoomed in) ────────────────────────────────────────

/// Walk the CIGAR, rendering actual bases at reference positions.
fn render_bases(
    read: &RenderRead,
    y: u16,
    area: Rect,
    transform: &ViewTransform,
    buf: &mut Buffer,
) {
    let mut read_pos: usize = 0;
    let mut ref_pos: u64 = read.start;

    for &op in &read.cigar_ops {
        match op {
            CigarOp::SoftClip(n) => {
                // Soft-clipped bases: not aligned to ref, advance read only
                read_pos += n as usize;
            }

            CigarOp::Match(n) | CigarOp::Mismatch(n) => {
                let is_mismatch_op = matches!(op, CigarOp::Mismatch(_));
                for _ in 0..n {
                    let base = read.sequence.get(read_pos).copied().unwrap_or(b'N');
                    let style = if is_mismatch_op {
                        // CIGAR explicitly marks mismatch
                        mismatch_style(base)
                    } else {
                        match_style(base, read.mapq, read.is_secondary || read.is_supplementary)
                    };
                    draw_at_ref_pos(ref_pos, y, base as char, style, area, transform, buf);
                    read_pos += 1;
                    ref_pos += 1;
                }
            }

            CigarOp::Insertion(n) => {
                // Mark the insertion with a caret on the left edge of where it falls
                let style = Style::default()
                    .fg(Color::Magenta)
                    .add_modifier(Modifier::BOLD);
                draw_at_ref_pos(ref_pos, y, '^', style, area, transform, buf);
                read_pos += n as usize;
                // insertion does not advance ref
            }

            CigarOp::Deletion(n) => {
                for _ in 0..n {
                    let style = Style::default().fg(Color::White).bg(Color::DarkGray);
                    draw_at_ref_pos(ref_pos, y, '-', style, area, transform, buf);
                    ref_pos += 1;
                }
            }

            CigarOp::Skip(n) => {
                // Intron / large skip: thin line
                for _ in 0..n {
                    let style = Style::default().fg(Color::DarkGray);
                    draw_at_ref_pos(ref_pos, y, '─', style, area, transform, buf);
                    ref_pos += 1;
                }
            }
        }
    }
}

/// Draw a single character at a reference position, translated to screen coordinates.
#[inline]
fn draw_at_ref_pos(
    ref_pos: u64,
    y: u16,
    ch: char,
    style: Style,
    area: Rect,
    transform: &ViewTransform,
    buf: &mut Buffer,
) {
    if ref_pos < transform.region_start || ref_pos >= transform.region_end {
        return;
    }
    let span = (transform.region_end - transform.region_start) as f64;
    let col = ((ref_pos - transform.region_start) as f64 / span * transform.cols as f64) as u16;
    let x = area.x + col;
    if x < area.x + area.width {
        if let Some(cell) = buf.cell_mut((x, y)) {
            cell.set_char(ch).set_style(style);
        }
    }
}

// ─── Arrow rendering (zoomed out) ────────────────────────────────────────────

fn render_arrows(
    read: &RenderRead,
    y: u16,
    area: Rect,
    transform: &ViewTransform,
    buf: &mut Buffer,
) {
    let (col_start, col_end) = transform.bp_range_to_cols(read.start, read.end);
    if col_start >= area.width || col_end == 0 {
        return;
    }
    let x_start = area.x + col_start;
    let x_end = (area.x + col_end).min(area.x + area.width);

    let mut ref_pos = read.start;
    for &op in &read.cigar_ops {
        let ref_len = op.ref_len();
        let op_end = ref_pos + ref_len.max(1);
        let (oc_start, oc_end) = transform.bp_range_to_cols(ref_pos, op_end);
        let ox_start = (area.x + oc_start).max(x_start);
        let ox_end = (area.x + oc_end).min(x_end);

        let (style, ch) = arrow_op_style(&op, read);
        for x in ox_start..ox_end {
            if x < area.x + area.width {
                if let Some(cell) = buf.cell_mut((x, y)) {
                    cell.set_char(ch).set_style(style);
                }
            }
        }
        if ref_len > 0 {
            ref_pos += ref_len;
        }
    }

    if read.cigar_ops.is_empty() && x_start < x_end {
        let style = mapq_style(read.mapq, read.is_secondary || read.is_supplementary);
        if let Some(cell) = buf.cell_mut((x_start, y)) {
            cell.set_char(strand_char(read.strand)).set_style(style);
        }
    }
}

fn arrow_op_style(op: &CigarOp, read: &RenderRead) -> (Style, char) {
    match op {
        CigarOp::Match(_) => (
            mapq_style(read.mapq, read.is_secondary || read.is_supplementary),
            strand_char(read.strand),
        ),
        CigarOp::Mismatch(_) => (
            Style::default().fg(Color::Red).add_modifier(Modifier::BOLD),
            'X',
        ),
        CigarOp::Insertion(_) => (Style::default().fg(Color::Magenta), '^'),
        CigarOp::Deletion(_) => (Style::default().fg(Color::White).bg(Color::DarkGray), '-'),
        CigarOp::Skip(_) => (Style::default().fg(Color::DarkGray), '─'),
        CigarOp::SoftClip(_) => (Style::default().fg(Color::DarkGray), '.'),
    }
}

// ─── Style helpers ────────────────────────────────────────────────────────────

/// Style for a matched base: color by nucleotide identity, dim if low mapq.
fn match_style(base: u8, mapq: u8, dim: bool) -> Style {
    let fg = base_color(base);
    // Background hint based on mapq so you can tell low-quality reads apart
    let bg = if mapq < 10 {
        Color::Reset
    } else {
        Color::Reset
    };
    let mut style = Style::default().fg(fg).bg(bg);
    if dim {
        style = style.add_modifier(Modifier::DIM);
    }
    style
}

/// Style for a CIGAR-X mismatch: bright, base-colored background.
fn mismatch_style(base: u8) -> Style {
    Style::default()
        .fg(Color::Black)
        .bg(base_color(base))
        .add_modifier(Modifier::BOLD)
}

fn mapq_style(mapq: u8, dim: bool) -> Style {
    let color = if mapq >= 60 {
        Color::White
    } else if mapq >= 30 {
        Color::Gray
    } else {
        Color::DarkGray
    };
    let mut s = Style::default().fg(color);
    if dim {
        s = s.add_modifier(Modifier::DIM);
    }
    s
}

/// Standard IGV-inspired nucleotide colors.
fn base_color(base: u8) -> Color {
    match base.to_ascii_uppercase() {
        b'A' => Color::Green,
        b'T' => Color::Red,
        b'G' => Color::Yellow,
        b'C' => Color::Blue,
        _ => Color::DarkGray,
    }
}

fn strand_char(strand: Strand) -> char {
    match strand {
        Strand::Forward => '>',
        Strand::Reverse => '<',
    }
}
