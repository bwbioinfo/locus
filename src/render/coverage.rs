use ratatui::{buffer::Buffer, layout::Rect, style::Style, widgets::Widget};

use crate::theme::Theme;

/// Unicode block characters for coverage histogram (8 levels + full block).
const BLOCKS: [char; 9] = [' ', '▁', '▂', '▃', '▄', '▅', '▆', '▇', '█'];

/// Coverage histogram widget.
pub struct CoverageTrack<'a> {
    pub bins: &'a [u32],
    pub theme: Theme,
}

impl<'a> Widget for CoverageTrack<'a> {
    fn render(self, area: Rect, buf: &mut Buffer) {
        if area.height == 0 || area.width == 0 || self.bins.is_empty() {
            return;
        }

        let max_cov = self.bins.iter().copied().max().unwrap_or(1).max(1);
        let rows = area.height as f64;

        for col in 0..area.width {
            // Map col to bin index
            let bin_idx = (col as usize * self.bins.len()) / area.width as usize;
            let bin_idx = bin_idx.min(self.bins.len() - 1);
            let cov = self.bins[bin_idx];

            // How many rows this column should fill
            let fill_frac = cov as f64 / max_cov as f64;
            let fill_rows = (fill_frac * rows).round() as u16;

            let col_style = coverage_style(fill_frac, self.theme);

            // Fill from the bottom
            for row_offset in 0..area.height {
                let row = area.height - 1 - row_offset;
                if row_offset < fill_rows {
                    if let Some(cell) = buf.cell_mut((area.x + col, area.y + row)) {
                        // Top partial block for the first filled cell
                        let ch = if row_offset + 1 == fill_rows && area.height > 1 {
                            let partial = ((fill_frac * rows).fract() * 8.0) as usize;
                            if partial > 0 {
                                BLOCKS[partial.min(8)]
                            } else {
                                BLOCKS[8]
                            }
                        } else {
                            BLOCKS[8]
                        };
                        cell.set_char(ch).set_style(col_style);
                    }
                }
            }
        }
    }
}

fn coverage_style(frac: f64, theme: Theme) -> Style {
    Style::default().fg(theme.coverage_color(frac))
}
