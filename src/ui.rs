use std::collections::BTreeMap;

use ratatui::{
    Frame,
    buffer::Buffer,
    layout::{Constraint, Direction, Layout, Rect},
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Clear, List, ListItem, Paragraph, Widget, Wrap},
};

use crate::render::{
    BASE_RENDER_THRESHOLD, ViewTransform,
    coverage::CoverageTrack,
    features::FeaturesTrack,
    reads::{ReadsTrack, SelectedPositionOverlay},
    reference::ReferenceTrack,
    ruler::Ruler,
};
use crate::{
    app::{App, Mode, ReadTrack},
    cache::{
        PhasePileupLayout, PhasePositionAlleleTallies, PileupRow, PositionAlleleTally, RenderRead,
    },
};

pub fn draw(frame: &mut Frame, app: &App) {
    let area = frame.area();
    let [top_bar, main, bottom_bar] = browser_layout(area);

    draw_top_bar(frame, app, top_bar);
    draw_main(frame, app, main);
    draw_bottom_bar(frame, app, bottom_bar);

    // Overlays (drawn on top)
    if app.show_help || app.mode == Mode::Help {
        draw_help_overlay(frame, app, area);
    }
    if app.mode == Mode::GoTo {
        draw_goto_overlay(frame, app, area);
    }
    if app.mode == Mode::FeatureSearch {
        draw_feature_search_overlay(frame, app, area);
    }
    if app.mode == Mode::ContigSelect {
        draw_contig_overlay(frame, app, area);
    }
    if app.mode == Mode::MapqFilter {
        draw_mapq_filter_overlay(frame, app, area);
    }
}

fn browser_layout(area: Rect) -> [Rect; 3] {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(1),
            Constraint::Min(5),
            Constraint::Length(1),
        ])
        .split(area);

    [chunks[0], chunks[1], chunks[2]]
}

const RULER_HEIGHT: u16 = 2;
const REFERENCE_HEIGHT: u16 = 1;
const FEATURES_HEIGHT: u16 = 4;
const MAX_COVERAGE_HEIGHT: u16 = 3;

fn coverage_height(main_height: u16) -> u16 {
    MAX_COVERAGE_HEIGHT.min(main_height / 5)
}

fn read_area_height(main_height: u16, has_reference: bool, has_features: bool) -> u16 {
    let fixed_height = RULER_HEIGHT
        .saturating_add(coverage_height(main_height))
        .saturating_add(if has_reference { REFERENCE_HEIGHT } else { 0 })
        .saturating_add(if has_features { FEATURES_HEIGHT } else { 0 });
    main_height.saturating_sub(fixed_height)
}

pub(crate) fn available_read_rows(
    terminal_rows: u16,
    has_reference: bool,
    has_features: bool,
) -> usize {
    let [_, main, _] = browser_layout(Rect::new(0, 0, 0, terminal_rows));
    read_area_height(main.height, has_reference, has_features) as usize
}

/// Return the genomic position under a terminal click in the main browser canvas.
pub(crate) fn genomic_position_at(app: &App, column: u16, row: u16) -> Option<u64> {
    let [_, main, _] = browser_layout(Rect::new(0, 0, app.terminal_cols, app.terminal_rows));
    if column < main.x
        || column >= main.x.saturating_add(main.width)
        || row < main.y
        || row >= main.y.saturating_add(main.height)
    {
        return None;
    }

    genomic_transform(app, main).col_to_bp(column.saturating_sub(main.x))
}

/// Return the read section under a terminal position for vertical navigation.
pub(crate) fn read_track_at(app: &App, column: u16, row: u16) -> Option<ReadTrack> {
    let [_, main, _] = browser_layout(Rect::new(0, 0, app.terminal_cols, app.terminal_rows));
    let reads_area = read_track_area(app, main);
    if !rect_contains(reads_area, column, row) {
        return None;
    }

    if !app.show_phasing {
        return Some(ReadTrack::Combined);
    }

    let layout = app.cache.phase_layout.as_ref()?;
    let [hp1, hp2, unphased] = phase_track_areas(reads_area, layout);
    if rect_contains(hp1, column, row) {
        Some(ReadTrack::Hp1)
    } else if rect_contains(hp2, column, row) {
        Some(ReadTrack::Hp2)
    } else if rect_contains(unphased, column, row) {
        Some(ReadTrack::Unphased)
    } else {
        None
    }
}

fn rect_contains(area: Rect, column: u16, row: u16) -> bool {
    column >= area.x
        && column < area.x.saturating_add(area.width)
        && row >= area.y
        && row < area.y.saturating_add(area.height)
}

fn draw_top_bar(frame: &mut Frame, app: &App, area: Rect) {
    let bp_per_col = app.view_span() as f64 / app.view_cols().max(1) as f64;
    let read_count = app.cache.pileup_rows.iter().map(Vec::len).sum::<usize>();
    let width = area.width as usize;
    let file_name = app
        .source
        .path
        .file_name()
        .unwrap_or_default()
        .to_string_lossy();

    let identity = format!(
        " LOCUS  file:{}  region:{}:{} ",
        file_name,
        app.current_contig(),
        format_region_display(app),
    );
    let insertion_mode = insertion_mode_label(app.expand_insertions);
    let selection_bracket_mode = selection_bracket_mode_label(app.show_selection_brackets);
    let methylation_mode = methylation_mode_label(app.show_methylation);
    let phasing_mode = phasing_mode_label(app.show_phasing);
    let theme_mode = theme_mode_label(app.theme);
    let mapq_filter = mapq_filter_label(app.min_mapq);
    let selected_info =
        selected_position_label(app.current_contig(), app.selected_ref_pos).map(|position| {
            let tally = if app.show_phasing {
                app.selected_phase_allele_tallies
                    .as_ref()
                    .map(phase_allele_tallies_label)
                    .unwrap_or_else(|| {
                        "HP1[none;m0/u0/r0] HP2[none;m0/u0/r0] U[none;m0/u0/r0]".to_string()
                    })
            } else {
                app.selected_allele_tally
                    .as_ref()
                    .map(selected_allele_tally_label)
                    .unwrap_or_else(|| "alleles:none meth:0 unmod:0 reads:0".to_string())
            };
            format!(" SEL {position}  {tally} ")
        });
    let metrics = format!(
        " {}  reads:{}  {}  scale:{:.1} bp/col  {}  {}  {}  {} ",
        mapq_filter,
        read_count,
        phasing_mode,
        bp_per_col,
        insertion_mode,
        selection_bracket_mode,
        methylation_mode,
        theme_mode
    );
    let status = app.status_msg.as_ref().map(|msg| format!(" status:{msg} "));
    let (identity, selected_info, metrics, status) = fit_top_bar(
        &identity,
        selected_info.as_deref(),
        &metrics,
        status.as_deref(),
        width,
    );

    let used = identity.len()
        + selected_info.as_ref().map_or(0, String::len)
        + metrics.len()
        + status.as_ref().map_or(0, String::len);
    let pad_len = width.saturating_sub(used);

    let mut spans = vec![Span::styled(
        identity,
        Style::default()
            .fg(app.theme.top_bar_identity_fg())
            .bg(app.theme.top_bar_identity_bg())
            .add_modifier(Modifier::BOLD),
    )];

    if let Some(selected_info) = selected_info {
        spans.push(Span::styled(
            selected_info,
            Style::default()
                .fg(app.theme.selected_info_fg())
                .bg(app.theme.selected_info_bg())
                .add_modifier(Modifier::BOLD),
        ));
    }

    spans.push(Span::raw(" ".repeat(pad_len)));
    spans.push(Span::styled(
        metrics,
        Style::default().fg(app.theme.top_bar_fg()),
    ));

    if let Some(status) = status {
        spans.push(Span::styled(
            status,
            Style::default().fg(app.theme.status_fg()),
        ));
    }

    frame.render_widget(
        Paragraph::new(Line::from(spans)).style(Style::default().bg(app.theme.top_bar_bg())),
        area,
    );
}

fn insertion_mode_label(expanded: bool) -> &'static str {
    if expanded {
        "ins:expanded"
    } else {
        "ins:collapsed"
    }
}

fn selection_bracket_mode_label(shown: bool) -> &'static str {
    if shown { "sel:brackets" } else { "sel:plain" }
}

fn methylation_mode_label(shown: bool) -> &'static str {
    if shown { "meth:on" } else { "meth:off" }
}

fn phasing_mode_label(shown: bool) -> &'static str {
    if shown { "phase:tracks" } else { "phase:off" }
}

fn theme_mode_label(theme: crate::theme::Theme) -> &'static str {
    match theme {
        crate::theme::Theme::Dark => "theme:dark",
        crate::theme::Theme::Light => "theme:light",
    }
}

fn mapq_filter_label(min_mapq: u8) -> String {
    if min_mapq == 0 {
        "mapq:all".to_string()
    } else {
        format!("mapq>={min_mapq}")
    }
}

fn selected_position_label(contig: &str, selected_ref_pos: Option<u64>) -> Option<String> {
    selected_ref_pos.map(|position| format!("{contig}:{}", position + 1))
}

fn selected_allele_tally_label(tally: &PositionAlleleTally) -> String {
    let alleles = allele_tally_details_label(tally);
    if alleles.is_empty() {
        format!(
            "alleles:none meth:{} unmod:{} reads:{}",
            tally.methylated_read_count, tally.unmodified_read_count, tally.total_read_count
        )
    } else {
        format!(
            "alleles:{alleles} meth:{} unmod:{} reads:{}",
            tally.methylated_read_count, tally.unmodified_read_count, tally.total_read_count
        )
    }
}

fn phase_allele_tallies_label(tallies: &PhasePositionAlleleTallies) -> String {
    format!(
        "HP1[{}] HP2[{}] U[{}]",
        phase_allele_tally_label(&tallies.hp1),
        phase_allele_tally_label(&tallies.hp2),
        phase_allele_tally_label(&tallies.unphased),
    )
}

fn phase_allele_tally_label(tally: &PositionAlleleTally) -> String {
    let alleles = compact_allele_tally_details_label(tally);
    if alleles.is_empty() {
        format!(
            "none;m{}/u{}/r{}",
            tally.methylated_read_count, tally.unmodified_read_count, tally.total_read_count
        )
    } else {
        format!(
            "{alleles};m{}/u{}/r{}",
            tally.methylated_read_count, tally.unmodified_read_count, tally.total_read_count
        )
    }
}

fn allele_tally_details_label(tally: &PositionAlleleTally) -> String {
    let mut alleles = Vec::new();
    for base in *b"ACGTN" {
        if let Some(count) = tally.base_counts.get(&base) {
            alleles.push(format!("{}:{count}", base as char));
        }
    }
    for (&base, &count) in &tally.base_counts {
        if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'N') {
            alleles.push(format!("{}:{count}", base as char));
        }
    }
    for (sequence, count) in &tally.deletion_counts {
        alleles.push(format!("-{}:{count}", String::from_utf8_lossy(sequence)));
    }
    if tally.deletion_count > 0 {
        alleles.push(format!("DEL:{}", tally.deletion_count));
    }
    for (sequence, count) in &tally.insertion_counts {
        alleles.push(format!("+{}:{count}", String::from_utf8_lossy(sequence)));
    }

    alleles.join(" ")
}

fn compact_allele_tally_details_label(tally: &PositionAlleleTally) -> String {
    let mut alleles = Vec::new();
    for base in *b"ACGTN" {
        if let Some(count) = tally.base_counts.get(&base) {
            alleles.push(format!("{}{}", base as char, count));
        }
    }
    for (&base, &count) in &tally.base_counts {
        if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'N') {
            alleles.push(format!("{}{}", base as char, count));
        }
    }
    for (sequence, count) in &tally.deletion_counts {
        alleles.push(format!("-{}{}", String::from_utf8_lossy(sequence), count));
    }
    if tally.deletion_count > 0 {
        alleles.push(format!("DEL{}", tally.deletion_count));
    }
    for (sequence, count) in &tally.insertion_counts {
        alleles.push(format!("+{}{}", String::from_utf8_lossy(sequence), count));
    }

    alleles.join(",")
}

fn truncate_to_width(text: &str, width: usize) -> String {
    if text.chars().count() <= width {
        return text.to_string();
    }
    if width == 0 {
        return String::new();
    }
    if width <= 1 {
        return "~".to_string();
    }

    let mut out = text
        .chars()
        .take(width.saturating_sub(1))
        .collect::<String>();
    out.push('~');
    out
}

fn fit_top_bar(
    identity: &str,
    selected_info: Option<&str>,
    metrics: &str,
    status: Option<&str>,
    width: usize,
) -> (String, Option<String>, String, Option<String>) {
    if width == 0 {
        return (
            String::new(),
            selected_info.map(|_| String::new()),
            String::new(),
            status.map(|_| String::new()),
        );
    }

    let identity_budget = if selected_info.is_some() && width >= 40 {
        width / 4
    } else if width < 40 {
        width / 2
    } else {
        width * 2 / 5
    };
    let identity = truncate_to_width(identity, identity_budget.max(1).min(width));
    let remaining = width.saturating_sub(identity.len());

    if let Some(selected_info) = selected_info {
        let selected_info = truncate_to_width(selected_info, remaining);
        let metrics = truncate_to_width(
            metrics,
            width.saturating_sub(identity.len() + selected_info.len()),
        );
        return (identity, Some(selected_info), metrics, None);
    }

    let status_reserve = status.map_or(0, |text| (remaining / 4).min(text.len()));
    let metrics = truncate_to_width(metrics, remaining.saturating_sub(status_reserve));
    let status = status
        .map(|text| truncate_to_width(text, width.saturating_sub(identity.len() + metrics.len())));

    (identity, None, metrics, status)
}

fn format_region_display(app: &App) -> String {
    format!("{}-{}", app.view_start + 1, app.view_end)
}

fn draw_main(frame: &mut Frame, app: &App, area: Rect) {
    let transform = genomic_transform(app, area);
    let chunks = main_chunks(app, area);
    let reference_h = if app.reference.is_some() {
        REFERENCE_HEIGHT
    } else {
        0
    };

    let mut chunk_idx = 0;

    // Ruler
    frame.render_widget(Ruler { transform }, chunks[chunk_idx]);
    chunk_idx += 1;

    if reference_h > 0 {
        frame.render_widget(
            ReferenceTrack {
                reference: app.cache.reference.as_ref(),
                transform,
                theme: app.theme,
            },
            chunks[chunk_idx],
        );
        chunk_idx += 1;
    }

    // Features track (only when GFF loaded)
    if let Some(ref gff) = app.gff {
        let visible = app.current_region();
        let feats = gff.features_in_region(&visible.contig, visible.start, visible.end);
        let feat_refs: Vec<&crate::gff::GffFeature> = feats.iter().collect();
        frame.render_widget(
            FeaturesTrack {
                features: &feat_refs,
                transform,
                theme: app.theme,
            },
            chunks[chunk_idx],
        );
        chunk_idx += 1;
    }

    // Coverage
    frame.render_widget(
        CoverageTrack {
            bins: &app.cache.coverage,
            theme: app.theme,
        },
        chunks[chunk_idx],
    );
    chunk_idx += 1;

    let reads_area = chunks[chunk_idx];
    if app.show_phasing {
        if let Some(layout) = app.cache.phase_layout.as_ref() {
            draw_phased_pileup(frame, app, transform, reads_area, layout);
        } else {
            draw_standard_pileup(frame, app, transform, reads_area);
        }
    } else {
        draw_standard_pileup(frame, app, transform, reads_area);
    }
}

fn main_chunks(app: &App, area: Rect) -> Vec<Rect> {
    let mut constraints = vec![Constraint::Length(RULER_HEIGHT)];
    if app.reference.is_some() {
        constraints.push(Constraint::Length(REFERENCE_HEIGHT));
    }
    if app.gff.is_some() {
        constraints.push(Constraint::Length(FEATURES_HEIGHT));
    }
    constraints.push(Constraint::Length(coverage_height(area.height)));
    constraints.push(Constraint::Min(read_area_height(
        area.height,
        app.reference.is_some(),
        app.gff.is_some(),
    )));

    Layout::default()
        .direction(Direction::Vertical)
        .constraints(constraints)
        .split(area)
        .to_vec()
}

fn read_track_area(app: &App, area: Rect) -> Rect {
    main_chunks(app, area).last().copied().unwrap_or_default()
}

fn genomic_transform(app: &App, area: Rect) -> ViewTransform {
    let base_transform =
        ViewTransform::new(app.view_start, app.view_end, area.width.saturating_sub(2));
    let insertion_gap = app.selected_insertion_gap(&base_transform);
    let show_selection_brackets =
        app.show_selection_brackets && base_transform.bp_per_col() <= BASE_RENDER_THRESHOLD;
    let selected_insertion = show_selection_brackets
        && app.selected_ref_pos.is_some_and(|selected_ref_pos| {
            insertion_gap.is_some_and(|gap| gap.anchor_ref_pos() == selected_ref_pos)
        });
    let selection_bracket = if show_selection_brackets && !selected_insertion {
        app.selected_ref_pos
    } else {
        None
    };
    base_transform
        .with_insertion_gap(insertion_gap)
        .with_selection_bracket(selection_bracket)
        .with_double_insertion_brackets(selected_insertion)
}

fn draw_standard_pileup(frame: &mut Frame, app: &App, transform: ViewTransform, area: Rect) {
    let offset = app
        .read_track_scroll(ReadTrack::Combined)
        .min(app.cache.pileup_rows.len());
    render_reads_track(
        frame,
        app,
        transform,
        &app.cache.pileup_rows[offset..],
        area,
    );
}

fn draw_phased_pileup(
    frame: &mut Frame,
    app: &App,
    transform: ViewTransform,
    area: Rect,
    layout: &PhasePileupLayout,
) {
    let areas = phase_track_areas(area, layout);
    draw_phase_section(
        frame,
        app,
        transform,
        PhaseSection {
            area: areas[0],
            track: ReadTrack::Hp1,
            label: "HP1",
            color: app.theme.phase_hp1_fg(),
            rows: &app.cache.pileup_rows[layout.hp1_rows.clone()],
            show_phase_set_boundaries: true,
        },
    );
    draw_phase_section(
        frame,
        app,
        transform,
        PhaseSection {
            area: areas[1],
            track: ReadTrack::Hp2,
            label: "HP2",
            color: app.theme.phase_hp2_fg(),
            rows: &app.cache.pileup_rows[layout.hp2_rows.clone()],
            show_phase_set_boundaries: true,
        },
    );

    if let Some(unphased_rows) = layout.unphased_rows.as_ref() {
        draw_phase_section(
            frame,
            app,
            transform,
            PhaseSection {
                area: areas[2],
                track: ReadTrack::Unphased,
                label: "Unphased",
                color: app.theme.phase_unphased_fg(),
                rows: &app.cache.pileup_rows[unphased_rows.clone()],
                show_phase_set_boundaries: false,
            },
        );
    }
}

struct PhaseSection<'a> {
    area: Rect,
    track: ReadTrack,
    label: &'a str,
    color: Color,
    rows: &'a [PileupRow],
    show_phase_set_boundaries: bool,
}

fn draw_phase_section(
    frame: &mut Frame,
    app: &App,
    transform: ViewTransform,
    section: PhaseSection<'_>,
) {
    let PhaseSection {
        area,
        track,
        label,
        color,
        rows,
        show_phase_set_boundaries,
    } = section;

    if area.height == 0 || area.width == 0 {
        return;
    }

    let phase_sets = if show_phase_set_boundaries {
        phase_set_boundaries(&app.cache.reads, rows)
    } else {
        Vec::new()
    };
    let read_count = rows.iter().map(Vec::len).sum::<usize>();
    let header = phase_section_header(label, read_count, &phase_sets, area.width as usize);
    let mut header_style = Style::default().fg(color).add_modifier(Modifier::BOLD);
    if app.active_read_track == track {
        header_style = header_style.add_modifier(Modifier::REVERSED);
    }
    frame.render_widget(
        Paragraph::new(header).style(header_style),
        Rect { height: 1, ..area },
    );

    render_phase_set_header_labels(
        frame,
        transform,
        &phase_sets,
        color,
        area,
        phase_section_prefix(label, read_count, &phase_sets)
            .chars()
            .count(),
    );

    let reads_area = Rect {
        y: area.y.saturating_add(1),
        height: area.height.saturating_sub(1),
        ..area
    };
    let offset = app.read_track_scroll(track).min(rows.len());
    render_reads_track(frame, app, transform, &rows[offset..], reads_area);
    frame.render_widget(
        PhaseSetBoundaryOverlay {
            boundaries: &phase_sets,
            transform,
            color,
        },
        reads_area,
    );
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PhaseSetBoundary {
    id: u32,
    start: u64,
}

fn phase_set_boundaries(reads: &[RenderRead], rows: &[PileupRow]) -> Vec<PhaseSetBoundary> {
    let mut starts = BTreeMap::new();

    for row in rows {
        for &read_idx in row {
            let Some(read) = reads.get(read_idx) else {
                continue;
            };
            let Some(id) = read.phase.phase_set else {
                continue;
            };
            starts
                .entry(id)
                .and_modify(|start: &mut u64| *start = (*start).min(read.start))
                .or_insert(read.start);
        }
    }

    let mut boundaries = starts
        .into_iter()
        .map(|(id, start)| PhaseSetBoundary { id, start })
        .collect::<Vec<_>>();
    boundaries.sort_by_key(|boundary| (boundary.start, boundary.id));
    boundaries
}

struct PhaseSetBoundaryOverlay<'a> {
    boundaries: &'a [PhaseSetBoundary],
    transform: ViewTransform,
    color: Color,
}

impl Widget for PhaseSetBoundaryOverlay<'_> {
    fn render(self, area: Rect, buf: &mut Buffer) {
        for boundary in self.boundaries {
            let Some(col) = self.transform.bp_to_col(boundary.start) else {
                continue;
            };
            let x = area.x.saturating_add(col);
            if x >= area.x.saturating_add(area.width) {
                continue;
            }

            for y in area.y..area.y.saturating_add(area.height) {
                let Some(cell) = buf.cell_mut((x, y)) else {
                    continue;
                };
                if cell.symbol().trim().is_empty() {
                    cell.set_symbol("┊");
                    cell.set_style(Style::default().fg(self.color).add_modifier(Modifier::DIM));
                } else {
                    cell.set_style(cell.style().add_modifier(Modifier::UNDERLINED));
                }
            }
        }
    }
}

fn render_phase_set_header_labels(
    frame: &mut Frame,
    transform: ViewTransform,
    boundaries: &[PhaseSetBoundary],
    color: Color,
    area: Rect,
    prefix_width: usize,
) {
    let mut previous_end = area.x.saturating_add(prefix_width as u16);
    let area_end = area.x.saturating_add(area.width);

    for boundary in boundaries {
        let Some(col) = transform.bp_to_col(boundary.start) else {
            continue;
        };
        let label = format!(" PS:{} ", boundary.id);
        let label_width = label.chars().count() as u16;
        let x = area.x.saturating_add(col);
        if x < previous_end || x.saturating_add(label_width) > area_end {
            continue;
        }

        frame.render_widget(
            Paragraph::new(label).style(Style::default().fg(color).add_modifier(Modifier::BOLD)),
            Rect {
                x,
                y: area.y,
                width: label_width,
                height: 1,
            },
        );
        previous_end = x.saturating_add(label_width);
    }
}

fn phase_section_prefix(label: &str, read_count: usize, phase_sets: &[PhaseSetBoundary]) -> String {
    let noun = if read_count == 1 { "read" } else { "reads" };
    let phase_sets = match phase_sets {
        [] => String::new(),
        [first] => format!("  PS:{}", first.id),
        [first, remaining @ ..] => format!("  PS:{} +{}", first.id, remaining.len()),
    };
    format!(" {label}  {read_count} {noun}{phase_sets} ")
}

fn phase_section_header(
    label: &str,
    read_count: usize,
    phase_sets: &[PhaseSetBoundary],
    width: usize,
) -> String {
    let prefix = phase_section_prefix(label, read_count, phase_sets);
    let divider = "─".repeat(width.saturating_sub(prefix.chars().count()));
    truncate_to_width(&format!("{prefix}{divider}"), width)
}

fn render_reads_track(
    frame: &mut Frame,
    app: &App,
    transform: ViewTransform,
    rows: &[PileupRow],
    area: Rect,
) {
    frame.render_widget(
        ReadsTrack {
            reads: &app.cache.reads,
            rows,
            reference: app.cache.reference.as_ref(),
            transform,
            show_names: area.width > 80,
            expand_insertions: app.expand_insertions,
            selected_ref_pos: app.selected_ref_pos,
            show_methylation: app.show_methylation,
            show_phasing: app.show_phasing,
            theme: app.theme,
        },
        area,
    );
    frame.render_widget(
        SelectedPositionOverlay {
            selected_ref_pos: app.selected_ref_pos,
            transform,
            theme: app.theme,
        },
        area,
    );
}

fn phase_track_areas(area: Rect, layout: &PhasePileupLayout) -> [Rect; 3] {
    let mut remaining_height = area.height;
    let hp1_height = phase_section_height(&mut remaining_height, layout.hp1_viewport_rows, true);
    let hp2_height = phase_section_height(&mut remaining_height, layout.hp2_viewport_rows, true);
    let unphased_height = phase_section_height(
        &mut remaining_height,
        layout.unphased_viewport_rows,
        layout.unphased_rows.is_some(),
    );

    let hp1 = Rect {
        height: hp1_height,
        ..area
    };
    let hp2 = Rect {
        y: area.y.saturating_add(hp1_height),
        height: hp2_height,
        ..area
    };
    let unphased = Rect {
        y: hp2.y.saturating_add(hp2_height),
        height: unphased_height,
        ..area
    };

    [hp1, hp2, unphased]
}

fn phase_section_height(remaining_height: &mut u16, row_count: usize, present: bool) -> u16 {
    if !present {
        return 0;
    }

    let requested = u16::try_from(row_count.saturating_add(1)).unwrap_or(u16::MAX);
    let height = requested.min(*remaining_height);
    *remaining_height = remaining_height.saturating_sub(height);
    height
}

fn draw_bottom_bar(frame: &mut Frame, app: &App, area: Rect) {
    let keys = match app.mode {
        Mode::Normal => {
            if app.gff.is_some() {
                " q:quit  ←/→:pan  Shift+←/→:1kb  Shift+↑/↓:scroll  Ctrl+↑/↓:track  +/-:zoom  i:insertions  b:brackets  m:methylation  p:phase tracks  Q:MAPQ  t:theme  Tab:next ins  g:goto  f:find  n/N:cycle  c:contigs  s:screenshot  Esc:clear  ?:help"
            } else {
                " q:quit  ←/→:pan  Shift+←/→:1kb  Shift+↑/↓:scroll  Ctrl+↑/↓:track  +/-:zoom  i:insertions  b:brackets  m:methylation  p:phase tracks  Q:MAPQ  t:theme  Tab:next ins  g:goto  c:contigs  r:refresh  s:screenshot  Esc:clear  ?:help"
            }
        }
        Mode::GoTo => " Enter:confirm  Esc:cancel",
        Mode::FeatureSearch => " type to search  Enter:jump  Tab/↑↓:cycle results  Esc:cancel",
        Mode::ContigSelect => " Enter:select  Esc:cancel",
        Mode::MapqFilter => " 0:show all  Enter:apply  Esc:cancel",
        Mode::Help => " Esc/q/?:close help",
    };
    frame.render_widget(
        Paragraph::new(keys).style(
            Style::default()
                .bg(app.theme.chrome_bg())
                .fg(app.theme.chrome_fg()),
        ),
        area,
    );
}

fn draw_goto_overlay(frame: &mut Frame, app: &App, area: Rect) {
    let popup = centered_rect(50, 12, area);
    let popup = Rect { height: 3, ..popup };
    frame.render_widget(Clear, popup);
    frame.render_widget(
        Paragraph::new(format!("Go to: {}_", app.command_buffer))
            .block(
                Block::default()
                    .title(" Jump to Region ")
                    .borders(Borders::ALL),
            )
            .style(Style::default().fg(app.theme.chrome_fg())),
        popup,
    );
}

fn draw_feature_search_overlay(frame: &mut Frame, app: &App, area: Rect) {
    let popup = centered_rect(60, 70, area);
    frame.render_widget(Clear, popup);

    // Split: input box on top, results list below
    let parts = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Length(3), Constraint::Min(3)])
        .split(popup);

    // Input
    let hint = if app.feature_matches.is_empty() && !app.command_buffer.is_empty() {
        " (no matches)"
    } else if !app.feature_matches.is_empty() {
        ""
    } else {
        " (type a gene or feature name)"
    };
    let title = format!(
        " Find Feature [{}/{}]{} ",
        if app.feature_matches.is_empty() {
            0
        } else {
            app.feature_match_idx + 1
        },
        app.feature_matches.len(),
        hint,
    );
    frame.render_widget(
        Paragraph::new(format!("{}_", app.command_buffer))
            .block(Block::default().title(title).borders(Borders::ALL))
            .style(Style::default().fg(app.theme.chrome_fg())),
        parts[0],
    );

    // Results list
    let gff = match app.gff.as_ref() {
        Some(g) => g,
        None => return,
    };

    let max_items = parts[1].height.saturating_sub(2) as usize;
    // Show a window around the current selection
    let total = app.feature_matches.len();
    let sel = app.feature_match_idx;
    let window_start = sel
        .saturating_sub(max_items / 2)
        .min(total.saturating_sub(max_items));

    let items: Vec<ListItem> = app
        .feature_matches
        .iter()
        .enumerate()
        .skip(window_start)
        .take(max_items)
        .map(|(i, &feat_idx)| {
            let feat = &gff.features[feat_idx];
            let marker = if i == sel { "▶ " } else { "  " };
            let style = if i == sel {
                Style::default()
                    .fg(app.theme.feature_label_fg())
                    .bg(app.theme.feature_color("gene"))
                    .add_modifier(Modifier::BOLD)
            } else {
                Style::default()
            };
            let coords = format!("{}:{}-{}", feat.seqname, feat.start + 1, feat.end);
            ListItem::new(format!(
                "{}{:<30} {:<12} {}",
                marker,
                feat.display_name(),
                feat.feature_type,
                coords,
            ))
            .style(style)
        })
        .collect();

    frame.render_widget(
        List::new(items).block(Block::default().title(" Results ").borders(Borders::ALL)),
        parts[1],
    );
}

fn draw_mapq_filter_overlay(frame: &mut Frame, app: &App, area: Rect) {
    let popup = centered_rect(50, 12, area);
    let popup = Rect { height: 3, ..popup };
    frame.render_widget(Clear, popup);
    frame.render_widget(
        Paragraph::new(format!(
            "Minimum MAPQ: {}_  (current: {})",
            app.command_buffer, app.min_mapq
        ))
        .block(
            Block::default()
                .title(" Read Quality Filter (0 shows all) ")
                .borders(Borders::ALL),
        )
        .style(Style::default().fg(app.theme.chrome_fg())),
        popup,
    );
}

fn draw_contig_overlay(frame: &mut Frame, app: &App, area: Rect) {
    let popup = centered_rect(40, 60, area);
    frame.render_widget(Clear, popup);

    let items: Vec<ListItem> = app
        .source
        .contigs
        .iter()
        .enumerate()
        .map(|(i, c)| {
            let marker = if i == app.contig_idx { "▶ " } else { "  " };
            let style = if i == app.contig_idx {
                Style::default()
                    .fg(app.theme.brand_fg())
                    .add_modifier(Modifier::BOLD)
            } else {
                Style::default()
            };
            ListItem::new(format!("{}{}: {} ({} bp)", marker, i + 1, c.name, c.length)).style(style)
        })
        .collect();

    frame.render_widget(
        List::new(items).block(
            Block::default()
                .title(" Select Contig (Enter number) ")
                .borders(Borders::ALL),
        ),
        popup,
    );
}

fn draw_help_overlay(frame: &mut Frame, app: &App, area: Rect) {
    let popup = centered_rect(60, 80, area);
    frame.render_widget(Clear, popup);

    let help_text = vec![
        Line::from(Span::styled(
            "  Locus Keybindings",
            Style::default().add_modifier(Modifier::BOLD),
        )),
        Line::from(""),
        Line::from("  q          Quit"),
        Line::from("  h / ←      Pan left (or move selected base left one bp)"),
        Line::from("  l / →      Pan right (or move selected base right one bp)"),
        Line::from("  Shift+←/→  Pan 1,000 bp left / right"),
        Line::from("  H / L      Pan 1,000 bp left / right"),
        Line::from("  ↑ / + / =  Zoom in"),
        Line::from("  ↓ / -      Zoom out"),
        Line::from("  Shift+↑/↓  Scroll active read track"),
        Line::from("  Ctrl+↑/↓   Select previous / next phased read track"),
        Line::from("  Mouse wheel Scroll read track under pointer"),
        Line::from("  Left click Select genomic position and highlight read bases"),
        Line::from("  Esc        Clear selected position"),
        Line::from("  b          Toggle fluorescent selection brackets"),
        Line::from("  i          Toggle expanded insertion sequence"),
        Line::from("  m          Toggle read methylation"),
        Line::from("  p          Toggle separated HP1 / HP2 read tracks"),
        Line::from("  Q          Set minimum read MAPQ (0 shows all)"),
        Line::from("  t          Toggle dark/light theme"),
        Line::from("  Tab        Move to next expanded insertion"),
        Line::from("  Shift+Tab  Move to previous expanded insertion"),
        Line::from("  g          Go to region  (e.g. chr1:1000-2000)"),
        Line::from("  f          Find feature / gene by name  (requires --gff)"),
        Line::from("  n / N      Cycle to next / previous feature match"),
        Line::from("  c          Contig selector"),
        Line::from("  r          Refresh current region"),
        Line::from("  s          Save ANSI text and HTML screenshots to screenshots/"),
        Line::from("  ?          Toggle this help"),
        Line::from(""),
        Line::from("  In feature search overlay:"),
        Line::from("    type     Filter results in real time"),
        Line::from("    Tab / ↓  Next result"),
        Line::from("    ↑        Previous result"),
        Line::from("    Enter    Jump to selected feature"),
        Line::from("    Esc      Close without jumping"),
        Line::from(""),
        Line::from("  Read colors:"),
        Line::from("    Phase off: MAPQ uses high / medium / low contrast"),
        Line::from("    Phase tracks: HP1 cyan, HP2 magenta, unphased gray"),
        Line::from("    MAPQ remains visible as bold / normal / dim intensity"),
        Line::from(
            "    Reference mismatches use base-colored bold backgrounds when --reference is loaded",
        ),
        Line::from("    The reversed phase-track header is the active keyboard target"),
        Line::from(""),
        Line::from("  CIGAR:  > / <  match   base highlight  mismatch   I  ins   -  del   ~  skip"),
        Line::from(""),
        Line::from("  Feature colors:"),
        Line::from("    Green  gene   Yellow  mRNA/transcript   Cyan  exon   Blue  CDS"),
        Line::from("    ─>─    intron/transcript backbone   █ exon   ▓ CDS   ▒ UTR"),
    ];

    frame.render_widget(
        Paragraph::new(help_text)
            .block(Block::default().title(" Help ").borders(Borders::ALL))
            .style(Style::default().fg(app.theme.chrome_fg()))
            .wrap(Wrap { trim: false }),
        popup,
    );
}

fn centered_rect(pct_x: u16, pct_y: u16, r: Rect) -> Rect {
    let w = r.width * pct_x / 100;
    let h = r.height * pct_y / 100;
    Rect {
        x: r.x + (r.width - w) / 2,
        y: r.y + (r.height - h) / 2,
        width: w,
        height: h,
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use ratatui::{Terminal, backend::TestBackend};

    use super::*;
    use crate::{
        bam::BamSource,
        cache::{CigarOp, ReadPhase, Strand},
        gff::GffStore,
        reference::ReferenceStore,
        region::Region,
        theme::Theme,
    };

    fn phase_read(name: &str, start: u64, phase_set: Option<u32>) -> RenderRead {
        RenderRead {
            name: name.to_string(),
            start,
            end: start + 10,
            strand: Strand::Forward,
            mapq: 60,
            cigar_ops: vec![CigarOp::Match(10)],
            sequence: b"AAAAAAAAAA".to_vec(),
            methylation: Vec::new(),
            deleted_reference_sequences: Vec::new(),
            phase: ReadPhase {
                haplotype: Some(1),
                phase_set,
            },
            is_secondary: false,
            is_supplementary: false,
            is_duplicate: false,
        }
    }

    #[test]
    fn truncate_to_width_respects_small_widths() {
        assert_eq!(truncate_to_width("abcdef", 0), "");
        assert_eq!(truncate_to_width("abcdef", 1), "~");
        assert_eq!(truncate_to_width("abcdef", 4), "abc~");
        assert_eq!(truncate_to_width("abc", 4), "abc");
        assert_eq!(truncate_to_width("──", 2), "──");
    }

    #[test]
    fn top_bar_preserves_identity_and_mapq_at_terminal_width() {
        let identity = " LOCUS  file:demo.sorted.bam  region:chrDemo:1-154 ";
        let metrics = " reads:3  mapq>=30  phase:tracks  scale:2.0 bp/col  ins:collapsed  meth:off  theme:dark ";
        let status = " status:minimum MAPQ set to 30 ";

        let (identity, selected_info, metrics, status) =
            fit_top_bar(identity, None, metrics, Some(status), 80);
        let status = status.expect("status remains present");

        assert!(identity.starts_with(" LOCUS"));
        assert!(selected_info.is_none());
        assert!(metrics.contains("mapq>=30"));
        assert!(metrics.contains("phase:tracks"));
        assert!(identity.len() + metrics.len() + status.len() <= 80);
    }

    #[test]
    fn top_bar_prioritizes_selected_information_and_mapq() {
        let identity = " LOCUS  file:demo.sorted.bam  region:chrDemo:45-115 ";
        let selected_info = " SEL chrDemo:60  alleles:C:2 T:2 +GGGG:1 ";
        let metrics = " mapq:all  reads:7  phase:tracks  scale:0.7 bp/col ";

        let (identity, selected_info, metrics, status) = fit_top_bar(
            identity,
            Some(selected_info),
            metrics,
            Some(" status:ignored "),
            110,
        );
        let selected_info = selected_info.expect("selection remains visible");

        assert!(identity.starts_with(" LOCUS"));
        assert!(selected_info.contains("SEL chrDemo:60"));
        assert!(selected_info.contains("+GGGG:1"));
        assert!(metrics.contains("mapq:all"));
        assert!(status.is_none());
        assert!(identity.len() + selected_info.len() + metrics.len() <= 110);
    }

    #[test]
    fn top_bar_keeps_all_phase_count_groups_at_demo_width() {
        let identity = " LOCUS  file:demo.sorted.bam  region:chrDemo:45-115 ";
        let selected_info = " SEL chrDemo:60  HP1[m2/u1/r3] HP2[m0/u0/r0] U[m1/u3/r4] ";
        let metrics = " mapq:all  reads:7  phase:tracks  scale:0.7 bp/col ";

        let (_, selected_info, _, _) =
            fit_top_bar(identity, Some(selected_info), metrics, None, 110);
        let selected_info = selected_info.expect("selection remains visible");

        assert!(selected_info.contains("HP1[m2/u1/r3]"));
        assert!(selected_info.contains("HP2[m0/u0/r0]"));
        assert!(selected_info.contains("U[m1/u3/r4]"));
    }

    #[test]
    fn methylation_mode_label_reflects_toggle_state() {
        assert_eq!(methylation_mode_label(false), "meth:off");
        assert_eq!(methylation_mode_label(true), "meth:on");
    }

    #[test]
    fn phasing_mode_label_reflects_toggle_state() {
        assert_eq!(phasing_mode_label(false), "phase:off");
        assert_eq!(phasing_mode_label(true), "phase:tracks");
    }

    #[test]
    fn phase_track_areas_are_compact_and_contiguous() {
        let area = Rect::new(7, 11, 80, 15);
        let layout = PhasePileupLayout {
            hp1_rows: 0..1,
            hp2_rows: 1..2,
            unphased_rows: Some(2..3),
            hp1_viewport_rows: 1,
            hp2_viewport_rows: 1,
            unphased_viewport_rows: 1,
        };
        let [hp1, hp2, unphased] = phase_track_areas(area, &layout);

        assert_eq!(hp1, Rect::new(7, 11, 80, 2));
        assert_eq!(hp2, Rect::new(7, 13, 80, 2));
        assert_eq!(unphased, Rect::new(7, 15, 80, 2));
        assert_eq!(unphased.y, hp2.y + hp2.height);
        assert!(unphased.y + unphased.height < area.y + area.height);
    }

    #[test]
    fn phase_track_areas_fit_dense_rows_and_omit_unphased_section() {
        let area = Rect::new(2, 3, 40, 7);
        let layout = PhasePileupLayout {
            hp1_rows: 0..3,
            hp2_rows: 3..5,
            hp1_viewport_rows: 3,
            hp2_viewport_rows: 2,
            ..PhasePileupLayout::default()
        };
        let [hp1, hp2, unphased] = phase_track_areas(area, &layout);

        assert_eq!(hp1, Rect::new(2, 3, 40, 4));
        assert_eq!(hp2, Rect::new(2, 7, 40, 3));
        assert_eq!(unphased, Rect::new(2, 10, 40, 0));
    }

    #[test]
    fn read_track_hit_testing_uses_the_phased_section_under_the_pointer() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/demo/demo.sorted.bam");
        let source = BamSource::open(path).expect("open demo BAM");
        let mut app = App::new(source, None, None, None, Theme::Dark, 0).expect("create app");
        app.terminal_cols = 80;
        app.terminal_rows = 24;
        app.show_phasing = true;
        app.cache.phase_layout = Some(PhasePileupLayout {
            hp1_rows: 0..3,
            hp2_rows: 3..6,
            unphased_rows: Some(6..8),
            hp1_viewport_rows: 2,
            hp2_viewport_rows: 2,
            unphased_viewport_rows: 1,
        });

        let [_, main, _] = browser_layout(Rect::new(0, 0, app.terminal_cols, app.terminal_rows));
        let reads_area = read_track_area(&app, main);
        let [hp1, hp2, unphased] = phase_track_areas(
            reads_area,
            app.cache.phase_layout.as_ref().expect("phase layout"),
        );

        assert_eq!(read_track_at(&app, hp1.x, hp1.y), Some(ReadTrack::Hp1));
        assert_eq!(read_track_at(&app, hp2.x, hp2.y), Some(ReadTrack::Hp2));
        assert_eq!(
            read_track_at(&app, unphased.x, unphased.y),
            Some(ReadTrack::Unphased)
        );
        assert_eq!(read_track_at(&app, 20, 1), None);
    }

    #[test]
    fn phase_section_height_reserves_a_header_for_empty_rows() {
        let mut remaining_height = 2;

        assert_eq!(phase_section_height(&mut remaining_height, 0, true), 1);
        assert_eq!(remaining_height, 1);
        assert_eq!(phase_section_height(&mut remaining_height, 0, false), 0);
    }

    #[test]
    fn phase_set_boundaries_are_sorted_deduplicated_and_ignore_missing_tags() {
        let reads = vec![
            phase_read("ps-100", 80, Some(100)),
            phase_read("ps-50", 50, Some(50)),
            phase_read("ps-100-earlier", 70, Some(100)),
            phase_read("untagged", 90, None),
        ];

        assert_eq!(
            phase_set_boundaries(&reads, &[vec![0, 1], vec![2, 3]]),
            vec![
                PhaseSetBoundary { id: 50, start: 50 },
                PhaseSetBoundary { id: 100, start: 70 },
            ]
        );
    }

    #[test]
    fn phase_set_boundary_overlay_marks_empty_cells_without_replacing_bases() {
        let area = Rect::new(0, 0, 4, 2);
        let transform = ViewTransform::new(100, 104, 4);
        let mut buffer = Buffer::empty(area);
        buffer[(1, 0)]
            .set_char('A')
            .set_style(Style::default().fg(Color::Green));

        PhaseSetBoundaryOverlay {
            boundaries: &[PhaseSetBoundary {
                id: 101,
                start: 101,
            }],
            transform,
            color: Color::Cyan,
        }
        .render(area, &mut buffer);

        assert_eq!(buffer[(1, 0)].symbol(), "A");
        assert_eq!(buffer[(1, 0)].style().fg, Some(Color::Green));
        assert!(
            buffer[(1, 0)]
                .style()
                .add_modifier
                .contains(Modifier::UNDERLINED)
        );
        assert_eq!(buffer[(1, 1)].symbol(), "┊");
        assert_eq!(buffer[(1, 1)].style().fg, Some(Color::Cyan));
    }

    #[test]
    fn rendered_phase_sections_follow_their_packed_rows() {
        let demo_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/demo");
        let source = BamSource::open(demo_dir.join("demo.sorted.bam")).expect("open demo BAM");
        let gff = GffStore::load(demo_dir.join("demo.sorted.gff.gz")).expect("open demo GFF");
        let reference = ReferenceStore::load(demo_dir.join("demo.fa")).expect("open demo FASTA");
        let mut app = App::new(
            source,
            Some(gff),
            Some(reference),
            Some(Region::new("chrDemo", 44, 115)),
            Theme::Dark,
            0,
        )
        .expect("create app");
        app.terminal_cols = 110;
        app.terminal_rows = 20;
        app.refresh().expect("load demo reads");
        app.toggle_phasing();

        let layout = app.cache.phase_layout.as_ref().expect("phased layout");
        assert_eq!(layout.hp1_rows.len(), 3);
        assert_eq!(layout.hp2_rows.len(), 1);
        assert_eq!(
            layout.unphased_rows.as_ref().map_or(0, |rows| rows.len()),
            1
        );

        let backend = TestBackend::new(app.terminal_cols, app.terminal_rows);
        let mut terminal = Terminal::new(backend).expect("test backend is infallible");
        terminal.draw(|frame| draw(frame, &app)).expect("draw app");

        let lines = (0..app.terminal_rows)
            .map(|row| {
                (0..app.terminal_cols)
                    .filter_map(|col| terminal.backend().buffer().cell((col, row)))
                    .map(|cell| cell.symbol())
                    .collect::<String>()
            })
            .collect::<Vec<_>>();
        let hp1 = lines
            .iter()
            .position(|line| line.contains("HP1  4 reads  PS:50 +1"))
            .expect("HP1 header");
        let hp2 = lines
            .iter()
            .position(|line| line.contains("HP2  1 read"))
            .expect("HP2 header");
        let unphased = lines
            .iter()
            .position(|line| line.contains("Unphased  2 reads"))
            .expect("unphased header");

        assert_eq!(hp2, hp1 + 4);
        assert_eq!(unphased, hp2 + 2);
        assert!(lines[hp1].contains("PS:100"));
        assert!(
            terminal
                .backend()
                .buffer()
                .cell((1, hp1 as u16))
                .expect("active header cell")
                .style()
                .add_modifier
                .contains(Modifier::REVERSED)
        );
    }

    #[test]
    fn phase_section_header_reports_group_and_phase_sets() {
        let phase_sets = [
            PhaseSetBoundary { id: 50, start: 50 },
            PhaseSetBoundary {
                id: 100,
                start: 100,
            },
        ];
        let header = phase_section_header("HP1", 4, &phase_sets, 40);

        assert!(header.starts_with(" HP1  4 reads  PS:50 +1 "));
        assert_eq!(header.chars().count(), 40);
        assert!(phase_section_header("HP2", 1, &[], 20).starts_with(" HP2  1 read "));
    }

    #[test]
    fn insertion_mode_label_reflects_toggle_state() {
        assert_eq!(insertion_mode_label(false), "ins:collapsed");
        assert_eq!(insertion_mode_label(true), "ins:expanded");
    }

    #[test]
    fn selection_bracket_mode_label_reflects_toggle_state() {
        assert_eq!(selection_bracket_mode_label(true), "sel:brackets");
        assert_eq!(selection_bracket_mode_label(false), "sel:plain");
    }

    #[test]
    fn theme_mode_label_reflects_theme_state() {
        assert_eq!(theme_mode_label(crate::theme::Theme::Dark), "theme:dark");
        assert_eq!(theme_mode_label(crate::theme::Theme::Light), "theme:light");
    }

    #[test]
    fn mapq_filter_label_reflects_threshold() {
        assert_eq!(mapq_filter_label(0), "mapq:all");
        assert_eq!(mapq_filter_label(30), "mapq>=30");
    }

    #[test]
    fn selected_position_label_is_one_based_and_includes_the_contig() {
        assert_eq!(
            selected_position_label("chrDemo", Some(64)),
            Some("chrDemo:65".to_string())
        );
        assert_eq!(selected_position_label("chrDemo", None), None);
    }

    #[test]
    fn selected_allele_tally_label_lists_bases_and_indels() {
        let tally = PositionAlleleTally {
            base_counts: BTreeMap::from([(b'A', 3), (b'C', 1)]),
            deletion_counts: BTreeMap::from([(b"ACT".to_vec(), 2)]),
            deletion_count: 1,
            insertion_counts: BTreeMap::from([(b"GG".to_vec(), 1)]),
            methylated_read_count: 2,
            unmodified_read_count: 3,
            total_read_count: 5,
        };

        assert_eq!(
            selected_allele_tally_label(&tally),
            "alleles:A:3 C:1 -ACT:2 DEL:1 +GG:1 meth:2 unmod:3 reads:5"
        );
        assert_eq!(
            selected_allele_tally_label(&PositionAlleleTally::default()),
            "alleles:none meth:0 unmod:0 reads:0"
        );
    }

    #[test]
    fn phase_allele_tallies_label_groups_all_selected_position_counts() {
        let tallies = PhasePositionAlleleTallies {
            hp1: PositionAlleleTally {
                base_counts: BTreeMap::from([(b'A', 3)]),
                methylated_read_count: 2,
                unmodified_read_count: 1,
                total_read_count: 3,
                ..PositionAlleleTally::default()
            },
            hp2: PositionAlleleTally {
                deletion_counts: BTreeMap::from([(b"G".to_vec(), 1)]),
                ..PositionAlleleTally::default()
            },
            unphased: PositionAlleleTally {
                insertion_counts: BTreeMap::from([(b"GG".to_vec(), 1)]),
                methylated_read_count: 1,
                unmodified_read_count: 3,
                total_read_count: 4,
                ..PositionAlleleTally::default()
            },
        };

        assert_eq!(
            phase_allele_tallies_label(&tallies),
            "HP1[A3;m2/u1/r3] HP2[-G1;m0/u0/r0] U[+GG1;m1/u3/r4]"
        );
    }

    #[test]
    fn browser_layout_reserves_the_top_and_bottom_bars() {
        let [top, main, bottom] = browser_layout(Rect::new(0, 0, 80, 24));

        assert_eq!(top, Rect::new(0, 0, 80, 1));
        assert_eq!(main, Rect::new(0, 1, 80, 22));
        assert_eq!(bottom, Rect::new(0, 23, 80, 1));
    }

    #[test]
    fn available_read_rows_matches_the_drawn_optional_tracks() {
        assert_eq!(available_read_rows(20, true, true), 8);
        assert_eq!(available_read_rows(20, false, false), 13);
    }

    #[test]
    fn insertion_gap_clicks_resolve_to_the_anchor_position() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/demo/demo.sorted.bam");
        let source = BamSource::open(path).expect("open demo BAM");
        let mut app = App::new(source, None, None, None, Theme::Dark, 0).expect("create app");
        app.terminal_cols = 110;
        app.terminal_rows = 20;
        app.jump_to_region(&Region::new("chrDemo", 44, 115))
            .expect("set demo region");
        app.refresh().expect("load demo reads");
        app.expand_insertions = true;
        app.cycle_insertion_expansion(true);
        let anchor = app.selected_insertion_ref_pos.expect("selected insertion");
        let [_, main, _] = browser_layout(Rect::new(0, 0, app.terminal_cols, app.terminal_rows));
        let transform = genomic_transform(&app, main);
        let gap = app
            .selected_insertion_gap(&transform)
            .expect("visible insertion gap");
        assert_eq!(anchor, gap.anchor_ref_pos());
        let (left_border, right_border) = transform
            .insertion_border_cols(gap.ref_pos)
            .expect("visible insertion gap");

        assert_eq!(
            genomic_position_at(&app, main.x + left_border, main.y),
            Some(anchor)
        );
        assert_eq!(
            genomic_position_at(&app, main.x + right_border, main.y),
            Some(anchor)
        );
    }

    #[test]
    fn selected_expanded_insertions_use_double_brackets_when_enabled() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/demo/demo.sorted.bam");
        let source = BamSource::open(path).expect("open demo BAM");
        let mut app = App::new(source, None, None, None, Theme::Dark, 0).expect("create app");
        app.terminal_cols = 110;
        app.terminal_rows = 20;
        app.jump_to_region(&Region::new("chrDemo", 44, 115))
            .expect("set demo region");
        app.refresh().expect("load demo reads");
        app.expand_insertions = true;
        app.cycle_insertion_expansion(true);
        let [_, main, _] = browser_layout(Rect::new(0, 0, app.terminal_cols, app.terminal_rows));
        let base_transform =
            ViewTransform::new(app.view_start, app.view_end, main.width.saturating_sub(2));
        let gap = app
            .selected_insertion_gap(&base_transform)
            .expect("visible insertion gap");

        app.select_reference_position(gap.anchor_ref_pos());
        let selected_transform = genomic_transform(&app, main);
        assert_eq!(selected_transform.insertion_bracket_count(), 2);
        assert_eq!(selected_transform.selection_bracket, None);

        app.toggle_selection_brackets();
        let unbracketed_transform = genomic_transform(&app, main);
        assert_eq!(unbracketed_transform.insertion_bracket_count(), 1);
    }

    #[test]
    fn selection_brackets_follow_the_base_render_threshold() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/demo/demo.sorted.bam");
        let source = BamSource::open(path).expect("open demo BAM");
        let mut app = App::new(source, None, None, None, Theme::Dark, 0).expect("create app");
        app.view_start = 0;
        app.selected_ref_pos = Some(2);

        app.view_end = 5;
        let base_visible = genomic_transform(&app, Rect::new(0, 0, 7, 1));
        assert_eq!(base_visible.bp_per_col(), BASE_RENDER_THRESHOLD);
        assert_eq!(base_visible.selection_bracket, Some(2));

        app.view_end = 6;
        let too_wide = genomic_transform(&app, Rect::new(0, 0, 7, 1));
        assert_eq!(too_wide.selection_bracket, None);

        app.view_end = 5;
        app.toggle_selection_brackets();
        let disabled = genomic_transform(&app, Rect::new(0, 0, 7, 1));
        assert_eq!(disabled.selection_bracket, None);
    }

    #[test]
    fn clicked_deletion_columns_select_the_deleted_reference_bases() {
        let demo_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/demo");
        let source = BamSource::open(demo_dir.join("demo.sorted.bam")).expect("open demo BAM");
        let reference = ReferenceStore::load(demo_dir.join("demo.fa")).expect("open demo FASTA");
        let mut app = App::new(
            source,
            None,
            Some(reference),
            Some(Region::new("chrDemo", 61, 65)),
            Theme::Dark,
            0,
        )
        .expect("create app");
        app.terminal_cols = 6;
        app.terminal_rows = 20;
        app.show_selection_brackets = false;
        app.refresh().expect("load demo reads");

        let [_, main, _] = browser_layout(Rect::new(0, 0, app.terminal_cols, app.terminal_rows));
        let transform = genomic_transform(&app, main);
        for deleted_position in 62..64 {
            let column = transform
                .bp_to_col(deleted_position)
                .expect("deletion column");
            let selected_position = genomic_position_at(&app, main.x + column, main.y)
                .expect("clickable deletion column");
            assert_eq!(selected_position, deleted_position);

            app.select_reference_position(selected_position);
            assert_eq!(
                app.selected_allele_tally
                    .as_ref()
                    .and_then(|tally| tally.deletion_counts.get(b"GT" as &[u8])),
                Some(&1)
            );
        }
    }
}
