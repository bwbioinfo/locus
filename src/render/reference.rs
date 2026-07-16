use ratatui::{
    buffer::Buffer,
    layout::Rect,
    style::{Modifier, Style},
    widgets::Widget,
};

use crate::{reference::ReferenceSlice, theme::Theme};

use super::ViewTransform;

const BASE_RENDER_THRESHOLD: f64 = 5.0;

pub struct ReferenceTrack<'a> {
    pub reference: Option<&'a ReferenceSlice>,
    pub transform: ViewTransform,
    pub theme: Theme,
}

impl Widget for ReferenceTrack<'_> {
    fn render(self, area: Rect, buf: &mut Buffer) {
        if area.height == 0 || area.width == 0 {
            return;
        }

        let y = area.y;
        let Some(reference) = self.reference else {
            render_label("ref unavailable", y, area, self.theme, buf);
            return;
        };

        if self.transform.bp_per_col() > BASE_RENDER_THRESHOLD {
            render_span(reference, y, area, &self.transform, self.theme, buf);
            return;
        }

        for ref_pos in self.transform.region_start..self.transform.region_end {
            let Some(base) = reference.base_at(ref_pos) else {
                continue;
            };
            let Some(col) = self.transform.bp_to_col(ref_pos) else {
                continue;
            };
            let x = area.x + col;
            if x >= area.x + area.width {
                continue;
            }
            if let Some(cell) = buf.cell_mut((x, y)) {
                cell.set_char(base as char)
                    .set_style(base_style(base, self.theme).add_modifier(Modifier::BOLD));
            }
        }
    }
}

fn render_span(
    reference: &ReferenceSlice,
    y: u16,
    area: Rect,
    transform: &ViewTransform,
    theme: Theme,
    buf: &mut Buffer,
) {
    let (col_start, col_end) = transform.bp_range_to_cols(reference.start, reference.end());
    let x_start = area.x + col_start;
    let x_end = (area.x + col_end).min(area.x + area.width);
    for x in x_start..x_end {
        if let Some(cell) = buf.cell_mut((x, y)) {
            cell.set_char('─')
                .set_style(Style::default().fg(theme.subtle_fg()));
        }
    }
}

fn render_label(label: &str, y: u16, area: Rect, theme: Theme, buf: &mut Buffer) {
    for (i, ch) in label.chars().take(area.width as usize).enumerate() {
        if let Some(cell) = buf.cell_mut((area.x + i as u16, y)) {
            cell.set_char(ch)
                .set_style(Style::default().fg(theme.subtle_fg()));
        }
    }
}

fn base_style(base: u8, theme: Theme) -> Style {
    Style::default().fg(theme.base_color(base))
}
