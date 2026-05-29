use ratatui::{
    buffer::Buffer,
    layout::Rect,
    style::{Color, Modifier, Style},
    widgets::Widget,
};

use crate::gff::GffFeature;

use super::ViewTransform;

/// Color scheme for feature types.
fn feature_color(ty: &str) -> Color {
    match ty {
        "gene" | "pseudogene" => Color::Green,
        "mRNA" | "transcript" | "lnc_RNA" | "ncRNA" | "pre_miRNA" => Color::Yellow,
        "exon" => Color::Cyan,
        "CDS" | "start_codon" | "stop_codon" => Color::Blue,
        "UTR" | "five_prime_UTR" | "three_prime_UTR" => Color::Magenta,
        _ => Color::DarkGray,
    }
}

/// Priority order for stacking rows: genes at top, then transcripts, then sub-features.
fn feature_priority(ty: &str) -> u8 {
    match ty {
        "gene" | "pseudogene" => 0,
        "mRNA" | "transcript" | "lnc_RNA" | "ncRNA" | "pre_miRNA" => 1,
        "exon" => 2,
        "CDS" => 3,
        _ => 4,
    }
}

pub struct FeaturesTrack<'a> {
    pub features: &'a [&'a GffFeature],
    pub transform: ViewTransform,
}

impl<'a> Widget for FeaturesTrack<'a> {
    fn render(self, area: Rect, buf: &mut Buffer) {
        if area.height == 0 || area.width == 0 || self.features.is_empty() {
            return;
        }

        // Lay features into rows (greedy packing, same algorithm as reads).
        let mut rows: Vec<Vec<&&GffFeature>> = Vec::new();
        let mut row_ends: Vec<u16> = Vec::new();

        // Sort by priority then start
        let mut sorted: Vec<&&GffFeature> = self.features.iter().collect();
        sorted.sort_by_key(|f| (feature_priority(f.feature_type.as_str()), f.start));

        for feat in &sorted {
            let (col_start, col_end) = self.transform.bp_range_to_cols(feat.start, feat.end);
            let target = row_ends
                .iter()
                .position(|&end| col_start >= end + 1)
                .unwrap_or(row_ends.len());

            if target >= area.height as usize {
                continue;
            }
            if target == rows.len() {
                rows.push(Vec::new());
                row_ends.push(0);
            }
            rows[target].push(feat);
            row_ends[target] = col_end;
        }

        for (row_idx, row) in rows.iter().enumerate() {
            let y = area.y + row_idx as u16;
            for feat in row {
                render_feature(feat, y, area, &self.transform, buf);
            }
        }
    }
}

fn render_feature(
    feat: &GffFeature,
    y: u16,
    area: Rect,
    transform: &ViewTransform,
    buf: &mut Buffer,
) {
    let (col_start, col_end) = transform.bp_range_to_cols(feat.start, feat.end);
    let x_start = area.x + col_start;
    let x_end = (area.x + col_end).min(area.x + area.width);
    if x_start >= x_end {
        return;
    }

    let color = feature_color(&feat.feature_type);
    let style = Style::default().fg(color);
    let width = (x_end - x_start) as usize;

    // Strand arrow char for the body
    use crate::cache::Strand;
    let body_ch = match feat.strand {
        Some(Strand::Forward) => '>',
        Some(Strand::Reverse) => '<',
        None => '─',
    };

    // Fill body
    for x in x_start..x_end {
        if let Some(cell) = buf.cell_mut((x, y)) {
            cell.set_char(body_ch).set_style(style);
        }
    }

    // Overlay name label if there's enough room (>= 3 cols)
    let name = feat.display_name();
    if width >= 3 && !name.is_empty() {
        let label: String = name.chars().take(width.saturating_sub(2)).collect();
        let label_x = x_start + 1;
        let label_style = Style::default()
            .fg(Color::Black)
            .bg(color)
            .add_modifier(Modifier::BOLD);
        for (i, ch) in label.chars().enumerate() {
            let lx = label_x + i as u16;
            if lx < x_end {
                if let Some(cell) = buf.cell_mut((lx, y)) {
                    cell.set_char(ch).set_style(label_style);
                }
            }
        }
    }
}
