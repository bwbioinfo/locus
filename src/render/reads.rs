use ratatui::{
    buffer::Buffer,
    layout::Rect,
    style::{Color, Modifier, Style},
    widgets::Widget,
};

use crate::cache::{CigarOp, PileupRow, RenderRead, Strand};
use crate::reference::ReferenceSlice;

use super::ViewTransform;

/// Show individual bases when the view is this many bp per column or narrower.
const BASE_RENDER_THRESHOLD: f64 = 5.0;
const INSERTION_EXPAND_THRESHOLD: f64 = 1.0;

pub struct ReadsTrack<'a> {
    pub reads: &'a [RenderRead],
    pub rows: &'a [PileupRow],
    pub reference: Option<&'a ReferenceSlice>,
    pub transform: ViewTransform,
    #[allow(dead_code)]
    pub show_names: bool,
    pub expand_insertions: bool,
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
                    render_bases(
                        read,
                        self.reference,
                        y,
                        area,
                        &self.transform,
                        self.expand_insertions,
                        buf,
                    );
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
    reference: Option<&ReferenceSlice>,
    y: u16,
    area: Rect,
    transform: &ViewTransform,
    expand_insertions: bool,
    buf: &mut Buffer,
) {
    if expand_insertions && transform.bp_per_col() <= INSERTION_EXPAND_THRESHOLD {
        render_bases_with_expanded_insertions(read, reference, y, area, transform, buf);
        return;
    }

    let mut read_pos: usize = 0;
    let mut ref_pos: u64 = read.start;
    let mut insertions = Vec::new();

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
                    let ref_base = reference.and_then(|reference| reference.base_at(ref_pos));
                    let style = aligned_base_style(
                        base,
                        ref_base,
                        is_mismatch_op,
                        read.mapq,
                        read.is_secondary || read.is_supplementary,
                    );
                    draw_at_ref_pos(ref_pos, y, base as char, style, area, transform, buf);
                    read_pos += 1;
                    ref_pos += 1;
                }
            }

            CigarOp::Insertion(n) => {
                insertions.push(ref_pos);
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

    for ref_pos in insertions {
        let style = insertion_style();
        draw_at_ref_pos(ref_pos, y, '^', style, area, transform, buf);
    }
}

/// Render base-level reads with inserted sequence expanded into extra screen columns.
///
/// The reference coordinate after an insertion is shifted to the right for this read row only.
/// This keeps inserted bases visible while preserving the reference-anchored left boundary.
fn render_bases_with_expanded_insertions(
    read: &RenderRead,
    reference: Option<&ReferenceSlice>,
    y: u16,
    area: Rect,
    transform: &ViewTransform,
    buf: &mut Buffer,
) {
    let mut read_pos: usize = 0;
    let mut ref_pos: u64 = read.start;
    let mut insertion_offset: u16 = 0;
    let mut pending_right_boundary = false;
    let dim = read.is_secondary || read.is_supplementary;

    for &op in &read.cigar_ops {
        match op {
            CigarOp::SoftClip(n) => {
                read_pos += n as usize;
            }
            CigarOp::Match(n) | CigarOp::Mismatch(n) => {
                let is_mismatch_op = matches!(op, CigarOp::Mismatch(_));
                for _ in 0..n {
                    let base = read.sequence.get(read_pos).copied().unwrap_or(b'N');
                    let ref_base = reference.and_then(|reference| reference.base_at(ref_pos));
                    let mut style =
                        aligned_base_style(base, ref_base, is_mismatch_op, read.mapq, dim);
                    if pending_right_boundary {
                        style = insertion_boundary_style(style);
                        pending_right_boundary = false;
                    }
                    draw_at_ref_pos_with_offset(
                        ref_pos,
                        insertion_offset,
                        y,
                        base as char,
                        style,
                        area,
                        transform,
                        buf,
                    );
                    read_pos += 1;
                    ref_pos += 1;
                }
            }
            CigarOp::Insertion(n) => {
                underline_left_insertion_boundary(
                    ref_pos,
                    insertion_offset,
                    y,
                    area,
                    transform,
                    buf,
                );

                for i in 0..n {
                    let base = read
                        .sequence
                        .get(read_pos + i as usize)
                        .copied()
                        .unwrap_or(b'N');
                    let style = insertion_base_style(base);
                    draw_at_ref_pos_with_offset(
                        ref_pos,
                        insertion_offset.saturating_add(i as u16),
                        y,
                        base as char,
                        style,
                        area,
                        transform,
                        buf,
                    );
                }

                read_pos += n as usize;
                insertion_offset = insertion_offset.saturating_add(n as u16);
                pending_right_boundary = true;
            }
            CigarOp::Deletion(n) => {
                for _ in 0..n {
                    let mut style = Style::default().fg(Color::White).bg(Color::DarkGray);
                    if pending_right_boundary {
                        style = insertion_boundary_style(style);
                        pending_right_boundary = false;
                    }
                    draw_at_ref_pos_with_offset(
                        ref_pos,
                        insertion_offset,
                        y,
                        '-',
                        style,
                        area,
                        transform,
                        buf,
                    );
                    ref_pos += 1;
                }
            }
            CigarOp::Skip(n) => {
                for _ in 0..n {
                    let mut style = Style::default().fg(Color::DarkGray);
                    if pending_right_boundary {
                        style = insertion_boundary_style(style);
                        pending_right_boundary = false;
                    }
                    draw_at_ref_pos_with_offset(
                        ref_pos,
                        insertion_offset,
                        y,
                        '─',
                        style,
                        area,
                        transform,
                        buf,
                    );
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

#[inline]
fn draw_at_ref_pos_with_offset(
    ref_pos: u64,
    col_offset: u16,
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
    let x = area.x + col.saturating_add(col_offset);
    if x < area.x + area.width {
        if let Some(cell) = buf.cell_mut((x, y)) {
            cell.set_char(ch).set_style(style);
        }
    }
}

fn underline_left_insertion_boundary(
    ref_pos: u64,
    insertion_offset: u16,
    y: u16,
    area: Rect,
    transform: &ViewTransform,
    buf: &mut Buffer,
) {
    let Some(left_ref_pos) = ref_pos.checked_sub(1) else {
        return;
    };
    if left_ref_pos < transform.region_start {
        return;
    }
    if let Some(cell) =
        cell_at_ref_pos_with_offset(left_ref_pos, insertion_offset, y, area, transform, buf)
    {
        let style = insertion_boundary_style(cell.style());
        cell.set_style(style);
    }
}

fn cell_at_ref_pos_with_offset<'a>(
    ref_pos: u64,
    col_offset: u16,
    y: u16,
    area: Rect,
    transform: &ViewTransform,
    buf: &'a mut Buffer,
) -> Option<&'a mut ratatui::buffer::Cell> {
    if ref_pos < transform.region_start || ref_pos >= transform.region_end {
        return None;
    }
    let span = (transform.region_end - transform.region_start) as f64;
    let col = ((ref_pos - transform.region_start) as f64 / span * transform.cols as f64) as u16;
    let x = area.x + col.saturating_add(col_offset);
    if x < area.x + area.width {
        buf.cell_mut((x, y))
    } else {
        None
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

fn insertion_style() -> Style {
    Style::default()
        .fg(Color::Magenta)
        .add_modifier(Modifier::BOLD)
}

fn insertion_base_style(base: u8) -> Style {
    Style::default()
        .fg(base_color(base))
        .bg(Color::DarkGray)
        .add_modifier(Modifier::BOLD)
}

fn insertion_boundary_style(style: Style) -> Style {
    style.add_modifier(Modifier::BOLD | Modifier::UNDERLINED)
}

/// Style for a matched base: color by nucleotide identity, dim if low mapq.
fn aligned_base_style(
    base: u8,
    ref_base: Option<u8>,
    is_mismatch_op: bool,
    mapq: u8,
    dim: bool,
) -> Style {
    if is_mismatch_op || ref_base.is_some_and(|ref_base| bases_mismatch(base, ref_base)) {
        mismatch_style(base)
    } else {
        match_style(base, mapq, dim)
    }
}

fn bases_mismatch(base: u8, ref_base: u8) -> bool {
    let base = base.to_ascii_uppercase();
    let ref_base = ref_base.to_ascii_uppercase();
    matches!(base, b'A' | b'C' | b'G' | b'T')
        && matches!(ref_base, b'A' | b'C' | b'G' | b'T')
        && base != ref_base
}

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bases_mismatch_only_for_canonical_base_differences() {
        assert!(!bases_mismatch(b'A', b'a'));
        assert!(bases_mismatch(b'A', b'C'));
        assert!(!bases_mismatch(b'N', b'A'));
        assert!(!bases_mismatch(b'A', b'N'));
    }

    #[test]
    fn reference_difference_uses_mismatch_style() {
        let style = aligned_base_style(b'A', Some(b'C'), false, 60, false);
        assert_eq!(style.bg, Some(Color::Green));

        let style = aligned_base_style(b'A', Some(b'A'), false, 60, false);
        assert_eq!(style.bg, Some(Color::Reset));
        assert_eq!(style.fg, Some(Color::Green));
    }

    #[test]
    fn expanded_insertions_render_sequence_between_underlined_boundaries() {
        let read = RenderRead {
            name: "read-with-ins".to_string(),
            start: 10,
            end: 13,
            strand: Strand::Forward,
            mapq: 60,
            cigar_ops: vec![CigarOp::Match(2), CigarOp::Insertion(2), CigarOp::Match(1)],
            sequence: b"ACGGG".to_vec(),
            is_secondary: false,
            is_supplementary: false,
            is_duplicate: false,
        };
        let area = Rect::new(0, 0, 8, 1);
        let transform = ViewTransform::new(10, 18, 8);
        let mut buf = Buffer::empty(area);

        render_bases(&read, None, 0, area, &transform, true, &mut buf);

        assert_eq!(buf[(0, 0)].symbol(), "A");
        assert_eq!(buf[(1, 0)].symbol(), "C");
        assert_eq!(buf[(2, 0)].symbol(), "G");
        assert_eq!(buf[(3, 0)].symbol(), "G");
        assert_eq!(buf[(4, 0)].symbol(), "G");
        assert!(
            buf[(1, 0)]
                .style()
                .add_modifier
                .contains(Modifier::BOLD | Modifier::UNDERLINED)
        );
        assert!(
            buf[(4, 0)]
                .style()
                .add_modifier
                .contains(Modifier::BOLD | Modifier::UNDERLINED)
        );
    }

    #[test]
    fn compact_insertions_are_overlaid_after_following_bases() {
        let read = RenderRead {
            name: "read-with-ins".to_string(),
            start: 10,
            end: 13,
            strand: Strand::Forward,
            mapq: 60,
            cigar_ops: vec![CigarOp::Match(2), CigarOp::Insertion(2), CigarOp::Match(1)],
            sequence: b"ACGGG".to_vec(),
            is_secondary: false,
            is_supplementary: false,
            is_duplicate: false,
        };
        let area = Rect::new(0, 0, 3, 1);
        let transform = ViewTransform::new(10, 25, 3);
        let mut buf = Buffer::empty(area);

        render_bases(&read, None, 0, area, &transform, false, &mut buf);

        assert_eq!(buf[(0, 0)].symbol(), "^");
    }

    #[test]
    fn expanded_insertions_keep_later_boundary_offsets() {
        let read = RenderRead {
            name: "read-with-two-ins".to_string(),
            start: 10,
            end: 14,
            strand: Strand::Forward,
            mapq: 60,
            cigar_ops: vec![
                CigarOp::Match(2),
                CigarOp::Insertion(1),
                CigarOp::Match(1),
                CigarOp::Insertion(1),
                CigarOp::Match(1),
            ],
            sequence: b"ACGTGA".to_vec(),
            is_secondary: false,
            is_supplementary: false,
            is_duplicate: false,
        };
        let area = Rect::new(0, 0, 10, 1);
        let transform = ViewTransform::new(10, 20, 10);
        let mut buf = Buffer::empty(area);

        render_bases(&read, None, 0, area, &transform, true, &mut buf);

        assert_eq!(buf[(0, 0)].symbol(), "A");
        assert_eq!(buf[(1, 0)].symbol(), "C");
        assert_eq!(buf[(2, 0)].symbol(), "G");
        assert_eq!(buf[(3, 0)].symbol(), "T");
        assert_eq!(buf[(4, 0)].symbol(), "G");
        assert_eq!(buf[(5, 0)].symbol(), "A");
        assert!(
            buf[(3, 0)]
                .style()
                .add_modifier
                .contains(Modifier::BOLD | Modifier::UNDERLINED)
        );
        assert!(
            buf[(5, 0)]
                .style()
                .add_modifier
                .contains(Modifier::BOLD | Modifier::UNDERLINED)
        );
    }

    #[test]
    fn insertions_are_collapsed_by_default_even_at_single_base_zoom() {
        let read = RenderRead {
            name: "read-with-ins".to_string(),
            start: 10,
            end: 13,
            strand: Strand::Forward,
            mapq: 60,
            cigar_ops: vec![CigarOp::Match(2), CigarOp::Insertion(2), CigarOp::Match(1)],
            sequence: b"ACGGG".to_vec(),
            is_secondary: false,
            is_supplementary: false,
            is_duplicate: false,
        };
        let area = Rect::new(0, 0, 8, 1);
        let transform = ViewTransform::new(10, 18, 8);
        let mut buf = Buffer::empty(area);

        render_bases(&read, None, 0, area, &transform, false, &mut buf);

        assert_eq!(buf[(0, 0)].symbol(), "A");
        assert_eq!(buf[(1, 0)].symbol(), "C");
        assert_eq!(buf[(2, 0)].symbol(), "^");
        assert_ne!(buf[(3, 0)].symbol(), "G");
    }
}
