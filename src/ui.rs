use ratatui::{
    Frame,
    layout::{Constraint, Direction, Layout, Rect},
    style::{Color, Modifier, Style},
    text::{Line, Span},
    widgets::{Block, Borders, Clear, List, ListItem, Paragraph, Wrap},
};

use crate::app::{App, Mode};
use crate::render::{
    ViewTransform, coverage::CoverageTrack, features::FeaturesTrack, reads::ReadsTrack,
    ruler::Ruler,
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
        draw_help_overlay(frame, area);
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
}

fn draw_top_bar(frame: &mut Frame, app: &App, area: Rect) {
    let bp_per_col = app.view_span() as f64 / app.view_cols().max(1) as f64;
    let read_count = app.cache.reads.len();
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
    let metrics = format!(" scale:{:.1} bp/col  reads:{} ", bp_per_col, read_count);
    let mut status = app.status_msg.as_ref().map(|msg| format!(" status:{msg} "));

    let status_width = status.as_ref().map_or(0, |s| s.len());
    let reserved = metrics.len() + status_width;
    let identity_width = width.saturating_sub(reserved);
    let identity = truncate_to_width(&identity, identity_width);
    let metrics = truncate_to_width(
        &metrics,
        width.saturating_sub(identity.len() + status_width),
    );

    if let Some(ref mut status) = status {
        let used = identity.len() + metrics.len();
        *status = truncate_to_width(status, width.saturating_sub(used));
    }

    let used = identity.len() + metrics.len() + status.as_ref().map_or(0, |s| s.len());
    let pad_len = width.saturating_sub(used);

    let mut spans = vec![
        Span::styled(
            identity,
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        ),
        Span::raw(" ".repeat(pad_len)),
        Span::styled(metrics, Style::default().fg(Color::White)),
    ];

    if let Some(status) = status {
        spans.push(Span::styled(status, Style::default().fg(Color::Yellow)));
    }

    frame.render_widget(
        Paragraph::new(Line::from(spans)).style(Style::default().bg(Color::DarkGray)),
        area,
    );
}

fn truncate_to_width(text: &str, width: usize) -> String {
    if text.len() <= width {
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

fn format_region_display(app: &App) -> String {
    format!("{}-{}", app.view_start + 1, app.view_end)
}

fn draw_main(frame: &mut Frame, app: &App, area: Rect) {
    let transform = ViewTransform::new(app.view_start, app.view_end, area.width.saturating_sub(2));

    let ruler_h = 2u16;
    let features_h: u16 = if app.gff.is_some() { 3 } else { 0 };
    let coverage_h = 3u16.min(area.height / 5);
    let reads_h = area
        .height
        .saturating_sub(ruler_h + features_h + coverage_h);

    let mut constraints = vec![Constraint::Length(ruler_h)];
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

    // Features track (only when GFF loaded)
    if let Some(ref gff) = app.gff {
        let visible = app.current_region();
        let feats: Vec<&crate::gff::GffFeature> = gff
            .features_in_region(&visible.contig, visible.start, visible.end)
            .collect();
        frame.render_widget(
            FeaturesTrack {
                features: &feats,
                transform,
            },
            chunks[chunk_idx],
        );
        chunk_idx += 1;
    }

    // Coverage
    frame.render_widget(
        CoverageTrack {
            bins: &app.cache.coverage,
        },
        chunks[chunk_idx],
    );
    chunk_idx += 1;

    // Reads pileup
    let show_names = area.width > 80;
    frame.render_widget(
        ReadsTrack {
            reads: &app.cache.reads,
            rows: &app.cache.pileup_rows,
            transform,
            show_names,
        },
        chunks[chunk_idx],
    );

    // Hidden reads notice
    if app.cache.hidden_reads > 0 {
        let msg = format!(" +{} reads hidden ", app.cache.hidden_reads);
        let notice_area = Rect {
            x: area.x,
            y: area.y + area.height.saturating_sub(1),
            width: msg.len() as u16,
            height: 1,
        };
        frame.render_widget(
            Paragraph::new(msg).style(Style::default().fg(Color::Yellow)),
            notice_area,
        );
    }
}

fn draw_bottom_bar(frame: &mut Frame, app: &App, area: Rect) {
    let keys = match app.mode {
        Mode::Normal => {
            if app.gff.is_some() {
                " q:quit  ←/→:pan  +/-:zoom  g:goto  f:find  n/N:cycle  c:contigs  s:screenshot  ?:help"
            } else {
                " q:quit  ←/→:pan  +/-:zoom  g:goto  c:contigs  r:refresh  s:screenshot  ?:help"
            }
        }
        Mode::GoTo => " Enter:confirm  Esc:cancel",
        Mode::FeatureSearch => " type to search  Enter:jump  Tab/↑↓:cycle results  Esc:cancel",
        Mode::ContigSelect => " Enter:select  Esc:cancel",
        Mode::Help => " Esc/q/?:close help",
    };
    frame.render_widget(
        Paragraph::new(keys).style(Style::default().bg(Color::DarkGray).fg(Color::White)),
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
            .style(Style::default().fg(Color::White)),
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
            .style(Style::default().fg(Color::White)),
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
                    .fg(Color::Black)
                    .bg(Color::Green)
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
                    .fg(Color::Cyan)
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

fn draw_help_overlay(frame: &mut Frame, area: Rect) {
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
        Line::from("    Green   MAPQ ≥ 60  Yellow  MAPQ ≥ 10  Gray  MAPQ < 10"),
        Line::from(""),
        Line::from("  CIGAR:  > / <  match   X  mismatch   I  ins   -  del   ~  skip   S  clip"),
        Line::from(""),
        Line::from("  Feature colors:"),
        Line::from("    Green  gene   Yellow  mRNA/transcript   Cyan  exon   Blue  CDS"),
    ];

    frame.render_widget(
        Paragraph::new(help_text)
            .block(Block::default().title(" Help ").borders(Borders::ALL))
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
    }
}
