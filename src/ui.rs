use ratatui::{
    Frame,
    layout::{Constraint, Direction, Layout, Rect},
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Clear, List, ListItem, Paragraph, Wrap},
};

use crate::render::{
    ViewTransform, coverage::CoverageTrack, features::FeaturesTrack, reads::ReadsTrack,
    reference::ReferenceTrack, ruler::Ruler,
};
use crate::{
    app::{App, Mode},
    cache::{PhasePileupLayout, PileupRow},
};

pub fn draw(frame: &mut Frame, app: &App) {
    let area = frame.area();

    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(1),
            Constraint::Min(5),
            Constraint::Length(1),
        ])
        .split(area);

    draw_top_bar(frame, app, chunks[0]);
    draw_main(frame, app, chunks[1]);
    draw_bottom_bar(frame, app, chunks[2]);

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
    let metrics = format!(
        " reads:{}  {}  {}  scale:{:.1} bp/col  {}  {}  {} ",
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
    let base_transform =
        ViewTransform::new(app.view_start, app.view_end, area.width.saturating_sub(2));
    let insertion_gap = app.selected_insertion_gap(&base_transform);
    let transform = base_transform.with_insertion_gap(insertion_gap);

    let ruler_h = 2u16;
    let reference_h: u16 = if app.reference.is_some() { 1 } else { 0 };
    let features_h: u16 = if app.gff.is_some() { 4 } else { 0 };
    let coverage_h = 3u16.min(area.height / 5);
    let reads_h = area
        .height
        .saturating_sub(ruler_h + reference_h + features_h + coverage_h);

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
    let areas = phase_track_areas(area, layout.unphased_rows.is_some());
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
    } = section;

    if area.height == 0 || area.width == 0 {
        return;
    }

    let read_count = rows.iter().map(Vec::len).sum::<usize>() + hidden_reads;
    let header = phase_section_header(label, read_count, hidden_reads, area.width as usize);
    frame.render_widget(
        Paragraph::new(header).style(Style::default().fg(color).add_modifier(Modifier::BOLD)),
        Rect { height: 1, ..area },
    );

    render_reads_track(
        frame,
        app,
        transform,
        rows,
        Rect {
            y: area.y.saturating_add(1),
            height: area.height.saturating_sub(1),
            ..area
        },
    );
}

fn phase_section_header(
    label: &str,
    read_count: usize,
    hidden_reads: usize,
    width: usize,
) -> String {
    let noun = if read_count == 1 { "read" } else { "reads" };
    let hidden = if hidden_reads > 0 {
        format!("  +{hidden_reads} hidden")
    } else {
        String::new()
    };
    let prefix = format!(" {label}  {read_count} {noun}{hidden} ");
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
}

fn phase_track_areas(area: Rect, has_unphased: bool) -> [Rect; 3] {
    let unphased_height = if has_unphased && area.height >= 4 {
        (area.height / 5).max(2)
    } else {
        0
    };
    let phased_height = area.height.saturating_sub(unphased_height);
    let hp1_height = phased_height.div_ceil(2);
    let hp2_height = phased_height / 2;

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
    use super::*;

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
    fn phase_track_areas_are_contiguous_and_bounded() {
        let area = Rect::new(7, 11, 80, 15);
        let [hp1, hp2, unphased] = phase_track_areas(area, true);

        assert_eq!(hp1, Rect::new(7, 11, 80, 6));
        assert_eq!(hp2, Rect::new(7, 17, 80, 6));
        assert_eq!(unphased, Rect::new(7, 23, 80, 3));
        assert_eq!(unphased.y + unphased.height, area.y + area.height);
    }

    #[test]
    fn phase_track_areas_split_only_haplotypes_without_unphased_reads() {
        let area = Rect::new(2, 3, 40, 7);
        let [hp1, hp2, unphased] = phase_track_areas(area, false);

        assert_eq!(hp1, Rect::new(2, 3, 40, 4));
        assert_eq!(hp2, Rect::new(2, 7, 40, 3));
        assert_eq!(unphased, Rect::new(2, 10, 40, 0));
    }

    #[test]
    fn phase_section_header_reports_group_and_hidden_reads() {
        let header = phase_section_header("HP1", 4, 2, 30);

        assert!(header.starts_with(" HP1  4 reads  +2 hidden "));
        assert_eq!(header.chars().count(), 30);
        assert!(phase_section_header("HP2", 1, 0, 20).starts_with(" HP2  1 read "));
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
}
