use ratatui::{
    buffer::Buffer,
    layout::Rect,
    style::{Color, Modifier, Style},
    widgets::Widget,
};

use crate::cache::{AlignedModifiedBaseCall, CigarOp, PileupRow, RenderRead, Strand};
use crate::reference::ReferenceSlice;
use crate::theme::Theme;

use super::{InsertionGap, ViewTransform};

/// Show individual bases when the view is this many bp per column or narrower.
const BASE_RENDER_THRESHOLD: f64 = 5.0;
const INSERTION_EXPAND_THRESHOLD: f64 = 1.0;

#[derive(Debug, Clone, Copy)]
struct InsertionEvent {
    ref_pos: u64,
    read_pos: usize,
    len: u64,
}

pub struct ReadsTrack<'a> {
    pub reads: &'a [RenderRead],
    pub rows: &'a [PileupRow],
    pub reference: Option<&'a ReferenceSlice>,
    pub transform: ViewTransform,
    #[allow(dead_code)]
    pub show_names: bool,
    pub expand_insertions: bool,
    pub show_methylation: bool,
    pub theme: Theme,
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
                    let context = BaseRenderContext {
                        reference: self.reference,
                        y,
                        area,
                        transform: &self.transform,
                        expand_insertions: self.expand_insertions,
                        show_methylation: self.show_methylation,
                        theme: self.theme,
                    };
                    render_bases(read, context, buf);
                } else {
                    render_arrows(read, y, area, &self.transform, self.theme, buf);
                }
            }
        }
    }
}

pub fn selected_insertion_gap(
    reads: &[RenderRead],
    rows: &[PileupRow],
    transform: &ViewTransform,
) -> Option<InsertionGap> {
    if transform.bp_per_col() > INSERTION_EXPAND_THRESHOLD {
        return None;
    }

    let gaps = visible_insertion_gaps(reads, rows, transform);
    let center = transform.region_start + (transform.region_end - transform.region_start) / 2;
    gaps.into_iter()
        .min_by_key(|gap| (gap.ref_pos.abs_diff(center), std::cmp::Reverse(gap.len)))
}

pub fn visible_insertion_gaps(
    reads: &[RenderRead],
    rows: &[PileupRow],
    transform: &ViewTransform,
) -> Vec<InsertionGap> {
    if transform.bp_per_col() > INSERTION_EXPAND_THRESHOLD {
        return Vec::new();
    }

    let mut gaps: Vec<InsertionGap> = Vec::new();
    for row in rows {
        for &read_idx in row {
            let Some(read) = reads.get(read_idx) else {
                continue;
            };
            for insertion in insertion_events(read) {
                if insertion.ref_pos < transform.region_start
                    || insertion.ref_pos >= transform.region_end
                {
                    continue;
                }
                match gaps.iter_mut().find(|gap| gap.ref_pos == insertion.ref_pos) {
                    Some(gap) => gap.len = gap.len.max(insertion.len),
                    None => gaps.push(InsertionGap {
                        ref_pos: insertion.ref_pos,
                        len: insertion.len,
                    }),
                }
            }
        }
    }

    gaps.sort_by_key(|gap| gap.ref_pos);
    gaps
}

// ─── Base-level rendering (zoomed in) ────────────────────────────────────────

#[derive(Clone, Copy)]
struct BaseRenderContext<'a> {
    reference: Option<&'a ReferenceSlice>,
    y: u16,
    area: Rect,
    transform: &'a ViewTransform,
    expand_insertions: bool,
    show_methylation: bool,
    theme: Theme,
}

/// Walk the CIGAR, rendering actual bases at reference positions.
fn render_bases(read: &RenderRead, context: BaseRenderContext<'_>, buf: &mut Buffer) {
    let mut read_pos: usize = 0;
    let mut ref_pos: u64 = read.start;
    let mut insertions = Vec::new();
    let dim = read.is_secondary || read.is_supplementary;
    let methylation_calls = if context.show_methylation {
        read.aligned_methylation()
    } else {
        Vec::new()
    };

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
                    let ref_base = context
                        .reference
                        .and_then(|reference| reference.base_at(ref_pos));
                    let mut style =
                        aligned_base_style(base, ref_base, is_mismatch_op, dim, context.theme);
                    if let Some(call) = methylation_call_at(&methylation_calls, read_pos, ref_pos) {
                        style = methylation_base_style(style, call.probability, context.theme);
                    }
                    draw_at_ref_pos(
                        ref_pos,
                        context.y,
                        base as char,
                        style,
                        context.area,
                        context.transform,
                        buf,
                    );
                    read_pos += 1;
                    ref_pos += 1;
                }
            }

            CigarOp::Insertion(n) => {
                insertions.push(InsertionEvent {
                    ref_pos,
                    read_pos,
                    len: n,
                });
                read_pos += n as usize;
                // insertion does not advance ref
            }

            CigarOp::Deletion(n) => {
                for _ in 0..n {
                    let style = deletion_style(context.theme);
                    draw_at_ref_pos(
                        ref_pos,
                        context.y,
                        '-',
                        style,
                        context.area,
                        context.transform,
                        buf,
                    );
                    ref_pos += 1;
                }
            }

            CigarOp::Skip(n) => {
                // Intron / large skip: thin line
                for _ in 0..n {
                    let style = skip_style(context.theme);
                    draw_at_ref_pos(
                        ref_pos,
                        context.y,
                        '─',
                        style,
                        context.area,
                        context.transform,
                        buf,
                    );
                    ref_pos += 1;
                }
            }
        }
    }

    let selected_gap = context
        .expand_insertions
        .then_some(context.transform.insertion_gap)
        .flatten();
    for insertion in insertions {
        if selected_gap.is_some_and(|gap| gap.ref_pos == insertion.ref_pos) {
            render_inserted_bases(read, insertion, context, buf);
        } else {
            draw_at_ref_pos(
                insertion.ref_pos,
                context.y,
                'I',
                insertion_marker_style(context.theme),
                context.area,
                context.transform,
                buf,
            );
        }
        emphasize_insertion_boundaries(
            insertion.ref_pos,
            context.y,
            context.area,
            context.transform,
            buf,
        );
    }

    if let Some(gap) = selected_gap
        && read.start < gap.ref_pos
        && read.end > gap.ref_pos
    {
        draw_insertion_box(gap.ref_pos, context, buf);
        emphasize_insertion_boundaries(
            gap.ref_pos,
            context.y,
            context.area,
            context.transform,
            buf,
        );
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
    let Some(col) = transform.bp_to_col(ref_pos) else {
        return;
    };
    let x = area.x + col;
    if x < area.x + area.width
        && let Some(cell) = buf.cell_mut((x, y))
    {
        cell.set_char(ch).set_style(style);
    }
}

fn render_inserted_bases(
    read: &RenderRead,
    insertion: InsertionEvent,
    context: BaseRenderContext<'_>,
    buf: &mut Buffer,
) {
    for i in 0..insertion.len {
        let base = read
            .sequence
            .get(insertion.read_pos + i as usize)
            .copied()
            .unwrap_or(b'N');
        let Some(col) = context.transform.insertion_col(insertion.ref_pos, i) else {
            continue;
        };
        let x = context.area.x + col;
        if x < context.area.x + context.area.width
            && let Some(cell) = buf.cell_mut((x, context.y))
        {
            cell.set_char(base as char)
                .set_style(insertion_base_style(base, context.theme));
        }
    }
}

fn draw_insertion_box(insertion_ref_pos: u64, context: BaseRenderContext<'_>, buf: &mut Buffer) {
    let Some((left_col, right_col)) = context.transform.insertion_border_cols(insertion_ref_pos)
    else {
        return;
    };
    let style = insertion_box_style(context.theme);
    let left_x = context.area.x + left_col;
    let right_x = context.area.x + right_col;

    if left_x < context.area.x + context.area.width
        && let Some(cell) = buf.cell_mut((left_x, context.y))
    {
        cell.set_char('[').set_style(style);
    }
    if right_x < context.area.x + context.area.width
        && let Some(cell) = buf.cell_mut((right_x, context.y))
    {
        cell.set_char(']').set_style(style);
    }
}

fn emphasize_insertion_boundaries(
    ref_pos: u64,
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
    if let Some(cell) = cell_at_ref_pos(left_ref_pos, y, area, transform, buf) {
        let style = insertion_boundary_style(cell.style());
        cell.set_style(style);
    }
    if let Some(cell) = cell_at_ref_pos(ref_pos, y, area, transform, buf) {
        let style = insertion_boundary_style(cell.style());
        cell.set_style(style);
    }
}

fn cell_at_ref_pos<'a>(
    ref_pos: u64,
    y: u16,
    area: Rect,
    transform: &ViewTransform,
    buf: &'a mut Buffer,
) -> Option<&'a mut ratatui::buffer::Cell> {
    let col = transform.bp_to_col(ref_pos)?;
    let x = area.x + col;
    if x < area.x + area.width {
        buf.cell_mut((x, y))
    } else {
        None
    }
}

fn insertion_events(read: &RenderRead) -> Vec<InsertionEvent> {
    let mut read_pos: usize = 0;
    let mut ref_pos: u64 = read.start;
    let mut insertions = Vec::new();

    for &op in &read.cigar_ops {
        match op {
            CigarOp::SoftClip(n) => read_pos += n as usize,
            CigarOp::Match(n) | CigarOp::Mismatch(n) => {
                read_pos += n as usize;
                ref_pos += n;
            }
            CigarOp::Insertion(n) => {
                insertions.push(InsertionEvent {
                    ref_pos,
                    read_pos,
                    len: n,
                });
                read_pos += n as usize;
            }
            CigarOp::Deletion(n) | CigarOp::Skip(n) => ref_pos += n,
        }
    }

    insertions
}

// ─── Arrow rendering (zoomed out) ────────────────────────────────────────────

fn render_arrows(
    read: &RenderRead,
    y: u16,
    area: Rect,
    transform: &ViewTransform,
    theme: Theme,
    buf: &mut Buffer,
) {
    let (col_start, col_end) = transform.bp_range_to_cols(read.start, read.end);
    if col_start >= area.width || col_end == 0 {
        return;
    }
    let x_start = area.x + col_start;
    let x_end = (area.x + col_end).min(area.x + area.width);

    let mut ref_pos = read.start;
    let mut insertions = Vec::new();
    for &op in &read.cigar_ops {
        if matches!(op, CigarOp::Insertion(_)) {
            insertions.push(ref_pos);
            continue;
        }

        let ref_len = op.ref_len();
        let op_end = ref_pos + ref_len.max(1);
        let (oc_start, oc_end) = transform.bp_range_to_cols(ref_pos, op_end);
        let ox_start = (area.x + oc_start).max(x_start);
        let ox_end = (area.x + oc_end).min(x_end);

        let (style, ch) = arrow_op_style(&op, read, theme);
        for x in ox_start..ox_end {
            if x < area.x + area.width
                && let Some(cell) = buf.cell_mut((x, y))
            {
                cell.set_char(ch).set_style(style);
            }
        }
        if ref_len > 0 {
            ref_pos += ref_len;
        }
    }

    for insertion_ref_pos in insertions {
        draw_at_ref_pos(
            insertion_ref_pos,
            y,
            'I',
            insertion_marker_style(theme),
            area,
            transform,
            buf,
        );
    }

    if read.cigar_ops.is_empty() && x_start < x_end {
        let style = mapq_style(read.mapq, read.is_secondary || read.is_supplementary, theme);
        if let Some(cell) = buf.cell_mut((x_start, y)) {
            cell.set_char(strand_char(read.strand)).set_style(style);
        }
    }
}

fn arrow_op_style(op: &CigarOp, read: &RenderRead, theme: Theme) -> (Style, char) {
    match op {
        CigarOp::Match(_) => (
            mapq_style(read.mapq, read.is_secondary || read.is_supplementary, theme),
            strand_char(read.strand),
        ),
        CigarOp::Mismatch(_) => (
            Style::default()
                .fg(theme.base_color(b'T'))
                .add_modifier(Modifier::BOLD),
            'X',
        ),
        CigarOp::Insertion(_) => (insertion_marker_style(theme), 'I'),
        CigarOp::Deletion(_) => (deletion_style(theme), '-'),
        CigarOp::Skip(_) => (skip_style(theme), '─'),
        CigarOp::SoftClip(_) => (skip_style(theme), '.'),
    }
}

// ─── Style helpers ────────────────────────────────────────────────────────────

fn insertion_marker_style(theme: Theme) -> Style {
    Style::default()
        .fg(theme.insertion_marker_fg())
        .bg(theme.insertion_marker_bg())
        .add_modifier(Modifier::BOLD | Modifier::REVERSED | Modifier::UNDERLINED)
}

fn insertion_base_style(base: u8, theme: Theme) -> Style {
    Style::default()
        .fg(base_color(base, theme))
        .bg(Color::Reset)
        .add_modifier(Modifier::BOLD)
}

fn insertion_box_style(theme: Theme) -> Style {
    Style::default()
        .fg(theme.insertion_marker_bg())
        .add_modifier(Modifier::BOLD)
}

fn insertion_boundary_style(style: Style) -> Style {
    style.add_modifier(Modifier::BOLD | Modifier::UNDERLINED)
}

fn methylation_call_at(
    calls: &[AlignedModifiedBaseCall],
    read_pos: usize,
    ref_pos: u64,
) -> Option<&crate::cache::ModifiedBaseCall> {
    calls
        .iter()
        .find(|call| call.call.read_pos == read_pos && call.ref_pos == Some(ref_pos))
        .map(|call| &call.call)
}

fn methylation_base_style(style: Style, probability: Option<u8>, theme: Theme) -> Style {
    match probability {
        Some(probability) if probability >= 192 => style
            .fg(theme.methylation_high_fg())
            .bg(theme.methylation_high_bg())
            .add_modifier(Modifier::BOLD | Modifier::UNDERLINED),
        Some(probability) if probability >= 128 => style
            .bg(theme.methylation_mid_bg())
            .add_modifier(Modifier::UNDERLINED),
        Some(_) => style
            .fg(theme.methylation_low_fg())
            .bg(theme.methylation_low_bg())
            .add_modifier(Modifier::UNDERLINED),
        None => style
            .fg(theme.methylation_high_bg())
            .add_modifier(Modifier::UNDERLINED),
    }
}

/// Style for a matched base: color by nucleotide identity, dim if low mapq.
fn aligned_base_style(
    base: u8,
    ref_base: Option<u8>,
    is_mismatch_op: bool,
    dim: bool,
    theme: Theme,
) -> Style {
    if is_mismatch_op || ref_base.is_some_and(|ref_base| bases_mismatch(base, ref_base)) {
        mismatch_style(base, theme)
    } else {
        match_style(base, dim, theme)
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
fn match_style(base: u8, dim: bool, theme: Theme) -> Style {
    let fg = base_color(base, theme);
    let mut style = Style::default().fg(fg).bg(Color::Reset);
    if dim {
        style = style.add_modifier(Modifier::DIM);
    }
    style
}

/// Style for a CIGAR-X mismatch: bright, base-colored background.
fn mismatch_style(base: u8, theme: Theme) -> Style {
    Style::default()
        .fg(theme.mismatch_fg())
        .bg(base_color(base, theme))
        .add_modifier(Modifier::BOLD)
}

fn mapq_style(mapq: u8, dim: bool, theme: Theme) -> Style {
    let color = if mapq >= 60 {
        theme.chrome_fg()
    } else if mapq >= 30 {
        theme.low_contrast_fg()
    } else {
        theme.subtle_fg()
    };
    let mut s = Style::default().fg(color);
    if dim {
        s = s.add_modifier(Modifier::DIM);
    }
    s
}

/// Standard IGV-inspired nucleotide colors.
fn deletion_style(theme: Theme) -> Style {
    Style::default().fg(theme.chrome_fg()).bg(theme.subtle_fg())
}

fn skip_style(theme: Theme) -> Style {
    Style::default().fg(theme.subtle_fg())
}

fn base_color(base: u8, theme: Theme) -> Color {
    theme.base_color(base)
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
    use crate::cache::{ModificationStrand, ModifiedBaseCall};

    fn render_test_bases(
        read: &RenderRead,
        area: Rect,
        transform: &ViewTransform,
        expand_insertions: bool,
        show_methylation: bool,
        theme: Theme,
        buf: &mut Buffer,
    ) {
        render_bases(
            read,
            BaseRenderContext {
                reference: None,
                y: 0,
                area,
                transform,
                expand_insertions,
                show_methylation,
                theme,
            },
            buf,
        );
    }

    #[test]
    fn bases_mismatch_only_for_canonical_base_differences() {
        assert!(!bases_mismatch(b'A', b'a'));
        assert!(bases_mismatch(b'A', b'C'));
        assert!(!bases_mismatch(b'N', b'A'));
        assert!(!bases_mismatch(b'A', b'N'));
    }

    #[test]
    fn reference_difference_uses_mismatch_style() {
        let style = aligned_base_style(b'A', Some(b'C'), false, false, Theme::Dark);
        assert_eq!(style.bg, Some(Color::Green));

        let style = aligned_base_style(b'A', Some(b'A'), false, false, Theme::Dark);
        assert_eq!(style.bg, Some(Color::Reset));
        assert_eq!(style.fg, Some(Color::Green));
    }

    #[test]
    fn light_theme_uses_readable_base_and_methylation_colors() {
        let base_style = aligned_base_style(b'G', Some(b'G'), false, false, Theme::Light);
        assert_eq!(base_style.fg, Some(Color::Rgb(132, 89, 0)));

        let methylated = methylation_base_style(base_style, Some(240), Theme::Light);
        assert_eq!(methylated.fg, Some(Color::White));
        assert_eq!(methylated.bg, Some(Color::Rgb(0, 102, 150)));
        assert!(
            methylated
                .add_modifier
                .contains(Modifier::BOLD | Modifier::UNDERLINED)
        );
    }

    #[test]
    fn expanded_insertions_open_shared_gap_at_selected_locus() {
        let read = RenderRead {
            name: "read-with-ins".to_string(),
            start: 10,
            end: 14,
            strand: Strand::Forward,
            mapq: 60,
            cigar_ops: vec![CigarOp::Match(2), CigarOp::Insertion(1), CigarOp::Match(2)],
            sequence: b"ACGTA".to_vec(),
            methylation: Vec::new(),
            is_secondary: false,
            is_supplementary: false,
            is_duplicate: false,
        };
        let area = Rect::new(0, 0, 8, 1);
        let transform = ViewTransform::new(10, 18, 8).with_insertion_gap(Some(InsertionGap {
            ref_pos: 12,
            len: 1,
        }));
        let mut buf = Buffer::empty(area);

        render_test_bases(&read, area, &transform, true, false, Theme::Dark, &mut buf);

        assert_eq!(buf[(0, 0)].symbol(), "A");
        assert_eq!(buf[(1, 0)].symbol(), "C");
        assert_eq!(buf[(2, 0)].symbol(), "[");
        assert_eq!(buf[(3, 0)].symbol(), "G");
        assert_eq!(buf[(4, 0)].symbol(), "]");
        assert_eq!(buf[(5, 0)].symbol(), "T");
        assert_eq!(buf[(6, 0)].symbol(), "A");
        assert!(
            buf[(1, 0)]
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
    fn compact_insertions_are_overlaid_after_following_bases() {
        let read = RenderRead {
            name: "read-with-ins".to_string(),
            start: 10,
            end: 13,
            strand: Strand::Forward,
            mapq: 60,
            cigar_ops: vec![CigarOp::Match(2), CigarOp::Insertion(2), CigarOp::Match(1)],
            sequence: b"ACGGG".to_vec(),
            methylation: Vec::new(),
            is_secondary: false,
            is_supplementary: false,
            is_duplicate: false,
        };
        let area = Rect::new(0, 0, 3, 1);
        let transform = ViewTransform::new(10, 25, 3);
        let mut buf = Buffer::empty(area);

        render_test_bases(&read, area, &transform, false, false, Theme::Dark, &mut buf);

        assert_eq!(buf[(0, 0)].symbol(), "I");
        assert_eq!(buf[(0, 0)].style().bg, Some(Color::Magenta));
    }

    #[test]
    fn visible_insertion_gaps_are_sorted_and_deduplicated() {
        let read_a = RenderRead {
            name: "read-a".to_string(),
            start: 10,
            end: 15,
            strand: Strand::Forward,
            mapq: 60,
            cigar_ops: vec![CigarOp::Match(4), CigarOp::Insertion(1), CigarOp::Match(1)],
            sequence: b"ACGTTA".to_vec(),
            methylation: Vec::new(),
            is_secondary: false,
            is_supplementary: false,
            is_duplicate: false,
        };
        let read_b = RenderRead {
            name: "read-b".to_string(),
            start: 10,
            end: 15,
            strand: Strand::Forward,
            mapq: 60,
            cigar_ops: vec![
                CigarOp::Match(2),
                CigarOp::Insertion(2),
                CigarOp::Match(2),
                CigarOp::Insertion(3),
                CigarOp::Match(1),
            ],
            sequence: b"ACGGTTAAAA".to_vec(),
            methylation: Vec::new(),
            is_secondary: false,
            is_supplementary: false,
            is_duplicate: false,
        };
        let reads = vec![read_a, read_b];
        let rows = vec![vec![0, 1]];
        let transform = ViewTransform::new(10, 20, 10);

        let gaps = visible_insertion_gaps(&reads, &rows, &transform);

        assert_eq!(
            gaps,
            vec![
                InsertionGap {
                    ref_pos: 12,
                    len: 2
                },
                InsertionGap {
                    ref_pos: 14,
                    len: 3
                },
            ]
        );
    }

    #[test]
    fn expanded_insertions_shift_by_selected_gap_only() {
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
            methylation: Vec::new(),
            is_secondary: false,
            is_supplementary: false,
            is_duplicate: false,
        };
        let area = Rect::new(0, 0, 10, 1);
        let transform = ViewTransform::new(10, 20, 10).with_insertion_gap(Some(InsertionGap {
            ref_pos: 12,
            len: 1,
        }));
        let mut buf = Buffer::empty(area);

        render_test_bases(&read, area, &transform, true, false, Theme::Dark, &mut buf);

        assert_eq!(buf[(0, 0)].symbol(), "A");
        assert_eq!(buf[(1, 0)].symbol(), "C");
        assert_eq!(buf[(2, 0)].symbol(), "[");
        assert_eq!(buf[(3, 0)].symbol(), "G");
        assert_eq!(buf[(4, 0)].symbol(), "]");
        assert_eq!(buf[(5, 0)].symbol(), "T");
        assert_eq!(buf[(6, 0)].symbol(), "I");
        assert_eq!(buf[(7, 0)].symbol(), " ");
        assert!(
            buf[(5, 0)]
                .style()
                .add_modifier
                .contains(Modifier::BOLD | Modifier::UNDERLINED)
        );
        assert!(
            buf[(6, 0)]
                .style()
                .add_modifier
                .contains(Modifier::BOLD | Modifier::UNDERLINED)
        );
    }

    #[test]
    fn shared_gap_shifts_reads_without_selected_insertion() {
        let read = RenderRead {
            name: "read-without-ins".to_string(),
            start: 10,
            end: 14,
            strand: Strand::Forward,
            mapq: 60,
            cigar_ops: vec![CigarOp::Match(4)],
            sequence: b"ACGT".to_vec(),
            methylation: Vec::new(),
            is_secondary: false,
            is_supplementary: false,
            is_duplicate: false,
        };
        let area = Rect::new(0, 0, 8, 1);
        let transform = ViewTransform::new(10, 18, 8).with_insertion_gap(Some(InsertionGap {
            ref_pos: 12,
            len: 2,
        }));
        let mut buf = Buffer::empty(area);

        render_test_bases(&read, area, &transform, true, false, Theme::Dark, &mut buf);

        assert_eq!(buf[(0, 0)].symbol(), "A");
        assert_eq!(buf[(1, 0)].symbol(), "C");
        assert_eq!(buf[(2, 0)].symbol(), "[");
        assert_eq!(buf[(3, 0)].symbol(), " ");
        assert_eq!(buf[(4, 0)].symbol(), " ");
        assert_eq!(buf[(5, 0)].symbol(), "]");
        assert_eq!(buf[(6, 0)].symbol(), "G");
        assert_eq!(buf[(7, 0)].symbol(), "T");
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
            methylation: Vec::new(),
            is_secondary: false,
            is_supplementary: false,
            is_duplicate: false,
        };
        let area = Rect::new(0, 0, 8, 1);
        let transform = ViewTransform::new(10, 18, 8);
        let mut buf = Buffer::empty(area);

        render_test_bases(&read, area, &transform, false, false, Theme::Dark, &mut buf);

        assert_eq!(buf[(0, 0)].symbol(), "A");
        assert_eq!(buf[(1, 0)].symbol(), "C");
        assert_eq!(buf[(2, 0)].symbol(), "I");
        assert_eq!(buf[(2, 0)].style().bg, Some(Color::Magenta));
        assert_ne!(buf[(3, 0)].symbol(), "G");
    }

    #[test]
    fn high_confidence_methylation_styles_aligned_base() {
        let read = RenderRead {
            name: "methylated-read".to_string(),
            start: 10,
            end: 14,
            strand: Strand::Forward,
            mapq: 60,
            cigar_ops: vec![CigarOp::Match(4)],
            sequence: b"ACGT".to_vec(),
            methylation: vec![ModifiedBaseCall {
                read_pos: 1,
                canonical_base: b'C',
                strand: ModificationStrand::Forward,
                modification: "m".to_string(),
                probability: Some(240),
            }],
            is_secondary: false,
            is_supplementary: false,
            is_duplicate: false,
        };
        let area = Rect::new(0, 0, 4, 1);
        let transform = ViewTransform::new(10, 14, 4);
        let mut buf = Buffer::empty(area);

        render_test_bases(&read, area, &transform, false, true, Theme::Dark, &mut buf);

        assert_eq!(buf[(1, 0)].symbol(), "C");
        assert_eq!(buf[(1, 0)].style().fg, Some(Color::Black));
        assert_eq!(buf[(1, 0)].style().bg, Some(Color::Cyan));
        assert!(
            buf[(1, 0)]
                .style()
                .add_modifier
                .contains(Modifier::BOLD | Modifier::UNDERLINED)
        );
    }

    #[test]
    fn low_confidence_methylation_styles_base_with_expanded_insertions() {
        let read = RenderRead {
            name: "low-methylated-read".to_string(),
            start: 10,
            end: 14,
            strand: Strand::Forward,
            mapq: 60,
            cigar_ops: vec![CigarOp::Match(2), CigarOp::Insertion(1), CigarOp::Match(2)],
            sequence: b"ACGTA".to_vec(),
            methylation: vec![ModifiedBaseCall {
                read_pos: 3,
                canonical_base: b'T',
                strand: ModificationStrand::Forward,
                modification: "m".to_string(),
                probability: Some(40),
            }],
            is_secondary: false,
            is_supplementary: false,
            is_duplicate: false,
        };
        let area = Rect::new(0, 0, 8, 1);
        let transform = ViewTransform::new(10, 18, 8).with_insertion_gap(Some(InsertionGap {
            ref_pos: 12,
            len: 1,
        }));
        let mut buf = Buffer::empty(area);

        render_test_bases(&read, area, &transform, true, true, Theme::Dark, &mut buf);

        assert_eq!(buf[(5, 0)].symbol(), "T");
        assert_eq!(buf[(5, 0)].style().fg, Some(Color::White));
        assert_eq!(buf[(5, 0)].style().bg, Some(Color::DarkGray));
        assert!(
            buf[(5, 0)]
                .style()
                .add_modifier
                .contains(Modifier::UNDERLINED)
        );
    }
}
