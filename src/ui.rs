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
    ViewTransform,
    coverage::CoverageTrack,
    features::FeaturesTrack,
    reads::{ReadsTrack, SelectedPositionOverlay},
    reference::ReferenceTrack,
    ruler::Ruler,
};
use crate::{
    app::{App, Mode},
    cache::{PhasePileupLayout, PileupRow, PositionAlleleTally, RenderRead},
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

fn draw_top_bar(frame: &mut Frame, app: &App, area: Rect) {
    let bp_per_col = app.view_span() as f64 / app.view_cols().max(1) as f64;
    let read_count =
        app.cache.pileup_rows.iter().map(Vec::len).sum::<usize>() + app.cache.hidden_reads;
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
    let methylation_mode = methylation_mode_label(app.show_methylation);
    let phasing_mode = phasing_mode_label(app.show_phasing);
    let theme_mode = theme_mode_label(app.theme);
    let mapq_filter = mapq_filter_label(app.min_mapq);
    let selected_position = selected_position_label(app.current_contig(), app.selected_ref_pos)
        .map(|position| {
            let tally = app
                .selected_allele_tally
                .as_ref()
                .map(selected_allele_tally_label)
                .unwrap_or_else(|| "alleles:none".to_string());
            format!(" pos:{position} {tally}")
        })
        .unwrap_or_default();
    let metrics = format!(
        "{selected_position}  reads:{}  {}  {}  scale:{:.1} bp/col  {}  {}  {} ",
        read_count,
        mapq_filter,
        phasing_mode,
        bp_per_col,
        insertion_mode,
        methylation_mode,
        theme_mode
    );
    let status = app.status_msg.as_ref().map(|msg| format!(" status:{msg} "));
    let (identity, metrics, status) = fit_top_bar(&identity, &metrics, status.as_deref(), width);

    let used = identity.len() + metrics.len() + status.as_ref().map_or(0, |s| s.len());
    let pad_len = width.saturating_sub(used);

    let mut spans = vec![
        Span::styled(
            identity,
            Style::default()
                .fg(app.theme.brand_fg())
                .add_modifier(Modifier::BOLD),
        ),
        Span::raw(" ".repeat(pad_len)),
        Span::styled(metrics, Style::default().fg(app.theme.chrome_fg())),
    ];

    if let Some(status) = status {
        spans.push(Span::styled(
            status,
            Style::default().fg(app.theme.status_fg()),
        ));
    }

    frame.render_widget(
        Paragraph::new(Line::from(spans)).style(Style::default().bg(app.theme.chrome_bg())),
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
    if tally.deletion_count > 0 {
        alleles.push(format!("DEL:{}", tally.deletion_count));
    }
    for (sequence, count) in &tally.insertion_counts {
        alleles.push(format!("+{}:{count}", String::from_utf8_lossy(sequence)));
    }

    if alleles.is_empty() {
        "alleles:none".to_string()
    } else {
        format!("alleles:{}", alleles.join(" "))
    }
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
    metrics: &str,
    status: Option<&str>,
    width: usize,
) -> (String, String, Option<String>) {
    if width == 0 {
        return (String::new(), String::new(), status.map(|_| String::new()));
    }

    let identity_budget = if width < 40 { width / 2 } else { width * 2 / 5 };
    let identity = truncate_to_width(identity, identity_budget.max(1).min(width));
    let remaining = width.saturating_sub(identity.len());
    let status_reserve = status.map_or(0, |text| (remaining / 4).min(text.len()));
    let metrics = truncate_to_width(metrics, remaining.saturating_sub(status_reserve));
    let status = status
        .map(|text| truncate_to_width(text, width.saturating_sub(identity.len() + metrics.len())));

    (identity, metrics, status)
}

fn format_region_display(app: &App) -> String {
    format!("{}-{}", app.view_start + 1, app.view_end)
}

fn draw_main(frame: &mut Frame, app: &App, area: Rect) {
    let transform = genomic_transform(app, area);

    let ruler_h = RULER_HEIGHT;
    let reference_h = if app.reference.is_some() {
        REFERENCE_HEIGHT
    } else {
        0
    };
    let features_h = if app.gff.is_some() {
        FEATURES_HEIGHT
    } else {
        0
    };
    let coverage_h = coverage_height(area.height);
    let reads_h = read_area_height(area.height, app.reference.is_some(), app.gff.is_some());

    let mut constraints = vec![Constraint::Length(ruler_h)];
    if reference_h > 0 {
        constraints.push(Constraint::Length(reference_h));
    }
    if features_h > 0 {
        constraints.push(Constraint::Length(features_h));
    }
    constraints.push(Constraint::Length(coverage_h));
    constraints.push(Constraint::Min(reads_h));

    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints(constraints)
        .split(area);

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

fn genomic_transform(app: &App, area: Rect) -> ViewTransform {
    let base_transform =
        ViewTransform::new(app.view_start, app.view_end, area.width.saturating_sub(2));
    base_transform.with_insertion_gap(app.selected_insertion_gap(&base_transform))
}

fn draw_standard_pileup(frame: &mut Frame, app: &App, transform: ViewTransform, area: Rect) {
    render_reads_track(frame, app, transform, &app.cache.pileup_rows, area);

    if app.cache.hidden_reads > 0 {
        let msg = format!(" +{} reads hidden ", app.cache.hidden_reads);
        let notice_area = Rect {
            x: area.x,
            y: area.y + area.height.saturating_sub(1),
            width: (msg.len() as u16).min(area.width),
            height: 1,
        };
        frame.render_widget(
            Paragraph::new(msg).style(Style::default().fg(app.theme.status_fg())),
            notice_area,
        );
    }
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
            label: "HP1",
            color: app.theme.phase_hp1_fg(),
            rows: &app.cache.pileup_rows[layout.hp1_rows.clone()],
            hidden_reads: layout.hp1_hidden,
            show_phase_set_boundaries: true,
        },
    );
    draw_phase_section(
        frame,
        app,
        transform,
        PhaseSection {
            area: areas[1],
            label: "HP2",
            color: app.theme.phase_hp2_fg(),
            rows: &app.cache.pileup_rows[layout.hp2_rows.clone()],
            hidden_reads: layout.hp2_hidden,
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
                label: "Unphased",
                color: app.theme.phase_unphased_fg(),
                rows: &app.cache.pileup_rows[unphased_rows.clone()],
                hidden_reads: layout.unphased_hidden,
                show_phase_set_boundaries: false,
            },
        );
    }
}

struct PhaseSection<'a> {
    area: Rect,
    label: &'a str,
    color: Color,
    rows: &'a [PileupRow],
    hidden_reads: usize,
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
        label,
        color,
        rows,
        hidden_reads,
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
    let read_count = rows.iter().map(Vec::len).sum::<usize>() + hidden_reads;
    let header = phase_section_header(
        label,
        read_count,
        hidden_reads,
        &phase_sets,
        area.width as usize,
    );
    frame.render_widget(
        Paragraph::new(header).style(Style::default().fg(color).add_modifier(Modifier::BOLD)),
        Rect { height: 1, ..area },
    );

    render_phase_set_header_labels(
        frame,
        transform,
        &phase_sets,
        color,
        area,
        phase_section_prefix(label, read_count, hidden_reads, &phase_sets)
            .chars()
            .count(),
    );

    let reads_area = Rect {
        y: area.y.saturating_add(1),
        height: area.height.saturating_sub(1),
        ..area
    };
    render_reads_track(frame, app, transform, rows, reads_area);
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

fn phase_section_prefix(
    label: &str,
    read_count: usize,
    hidden_reads: usize,
    phase_sets: &[PhaseSetBoundary],
) -> String {
    let noun = if read_count == 1 { "read" } else { "reads" };
    let hidden = if hidden_reads > 0 {
        format!("  +{hidden_reads} hidden")
    } else {
        String::new()
    };
    let phase_sets = match phase_sets {
        [] => String::new(),
        [first] => format!("  PS:{}", first.id),
        [first, remaining @ ..] => format!("  PS:{} +{}", first.id, remaining.len()),
    };
    format!(" {label}  {read_count} {noun}{hidden}{phase_sets} ")
}

fn phase_section_header(
    label: &str,
    read_count: usize,
    hidden_reads: usize,
    phase_sets: &[PhaseSetBoundary],
    width: usize,
) -> String {
    let prefix = phase_section_prefix(label, read_count, hidden_reads, phase_sets);
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
        },
        area,
    );
}

fn phase_track_areas(area: Rect, layout: &PhasePileupLayout) -> [Rect; 3] {
    let mut remaining_height = area.height;
    let hp1_height = phase_section_height(&mut remaining_height, layout.hp1_rows.len(), true);
    let hp2_height = phase_section_height(&mut remaining_height, layout.hp2_rows.len(), true);
    let unphased_height = phase_section_height(
        &mut remaining_height,
        layout.unphased_rows.as_ref().map_or(0, |rows| rows.len()),
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
                " q:quit  ←/→:pan  +/-:zoom  i:insertions  m:methylation  p:phase tracks  Q:MAPQ  t:theme  Tab:next ins  g:goto  f:find  n/N:cycle  c:contigs  s:screenshot  ?:help"
            } else {
                " q:quit  ←/→:pan  +/-:zoom  i:insertions  m:methylation  p:phase tracks  Q:MAPQ  t:theme  Tab:next ins  g:goto  c:contigs  r:refresh  s:screenshot  ?:help"
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
        Line::from("  h / ←      Pan left (small)"),
        Line::from("  l / →      Pan right (small)"),
        Line::from("  H          Pan left (large)"),
        Line::from("  L          Pan right (large)"),
        Line::from("  ↑ / + / =  Zoom in"),
        Line::from("  ↓ / -      Zoom out"),
        Line::from("  Left click Select genomic position and highlight read bases"),
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

        let (identity, metrics, status) = fit_top_bar(identity, metrics, Some(status), 80);
        let status = status.expect("status remains present");

        assert!(identity.starts_with(" LOCUS"));
        assert!(metrics.contains("mapq>=30"));
        assert!(metrics.contains("phase:tracks"));
        assert!(identity.len() + metrics.len() + status.len() <= 80);
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
            ..PhasePileupLayout::default()
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
            ..PhasePileupLayout::default()
        };
        let [hp1, hp2, unphased] = phase_track_areas(area, &layout);

        assert_eq!(hp1, Rect::new(2, 3, 40, 4));
        assert_eq!(hp2, Rect::new(2, 7, 40, 3));
        assert_eq!(unphased, Rect::new(2, 10, 40, 0));
    }

    #[test]
    fn phase_section_height_reserves_a_header_for_hidden_rows() {
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
    }

    #[test]
    fn phase_section_header_reports_group_and_hidden_reads() {
        let phase_sets = [
            PhaseSetBoundary { id: 50, start: 50 },
            PhaseSetBoundary {
                id: 100,
                start: 100,
            },
        ];
        let header = phase_section_header("HP1", 4, 2, &phase_sets, 40);

        assert!(header.starts_with(" HP1  4 reads  +2 hidden  PS:50 +1 "));
        assert_eq!(header.chars().count(), 40);
        assert!(phase_section_header("HP2", 1, 0, &[], 20).starts_with(" HP2  1 read "));
    }

    #[test]
    fn insertion_mode_label_reflects_toggle_state() {
        assert_eq!(insertion_mode_label(false), "ins:collapsed");
        assert_eq!(insertion_mode_label(true), "ins:expanded");
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
            deletion_count: 2,
            insertion_counts: BTreeMap::from([(b"GG".to_vec(), 1)]),
        };

        assert_eq!(
            selected_allele_tally_label(&tally),
            "alleles:A:3 C:1 DEL:2 +GG:1"
        );
        assert_eq!(
            selected_allele_tally_label(&PositionAlleleTally::default()),
            "alleles:none"
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
        let (left_border, right_border) = transform
            .insertion_border_cols(anchor)
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
}
