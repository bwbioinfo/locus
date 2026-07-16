use ratatui::{
    buffer::Buffer,
    layout::Rect,
    style::{Color, Style},
    widgets::Widget,
};

use super::ViewTransform;

/// Coordinate ruler widget.
pub struct Ruler {
    pub transform: ViewTransform,
}

impl Widget for Ruler {
    fn render(self, area: Rect, buf: &mut Buffer) {
        if area.height == 0 || area.width == 0 {
            return;
        }

        let span = (self.transform.region_end - self.transform.region_start) as f64;
        let cols = area.width as f64;
        let bp_per_col = span / cols;

        // Pick a sensible tick interval (round to a power-of-10 or half)
        let tick_every = nice_interval(bp_per_col * 10.0);

        let tick_style = Style::default().fg(Color::DarkGray);
        let label_style = Style::default().fg(Color::White);

        // Ruler line (top row of the area)
        let y = area.y;
        for col in 0..area.width {
            let bp = self.transform.region_start + (col as f64 * bp_per_col) as u64;
            let cell = buf.cell_mut((area.x + col, y));
            if let Some(cell) = cell {
                // Display 1-based coordinate
                let pos_1based = bp + 1;
                if tick_every > 0 && pos_1based.is_multiple_of(tick_every) {
                    cell.set_char('|').set_style(tick_style);
                } else {
                    cell.set_char('─').set_style(tick_style);
                }
            }
        }

        // Label row (second row if available)
        if area.height >= 2 {
            let label_y = area.y + 1;
            let mut last_label_end = 0u16;

            for col in 0..area.width {
                let bp = self.transform.region_start + (col as f64 * bp_per_col) as u64;
                let pos_1based = bp + 1;
                if tick_every > 0 && pos_1based.is_multiple_of(tick_every) {
                    let label = format_position(pos_1based);
                    // Only print if it won't overlap the previous label
                    if col >= last_label_end {
                        for (i, ch) in label.chars().enumerate() {
                            let lc = col + i as u16;
                            if lc < area.width
                                && let Some(cell) = buf.cell_mut((area.x + lc, label_y))
                            {
                                cell.set_char(ch).set_style(label_style);
                            }
                        }
                        last_label_end = col + label.len() as u16 + 1;
                    }
                }
            }
        }
    }
}

/// Format a genomic position with comma separators, e.g. 1_234_567 → "1,234,567"
fn format_position(pos: u64) -> String {
    let s = pos.to_string();
    let mut out = String::new();
    for (i, ch) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            out.push(',');
        }
        out.push(ch);
    }
    out.chars().rev().collect()
}

/// Pick a "nice" tick interval >= `min_bp` (power of 10 or 5x/2x multiples).
fn nice_interval(min_bp: f64) -> u64 {
    if min_bp <= 0.0 {
        return 1;
    }
    let magnitude = 10f64.powf(min_bp.log10().floor());
    for &factor in &[1.0f64, 2.0, 5.0, 10.0] {
        let candidate = (magnitude * factor) as u64;
        if candidate as f64 >= min_bp {
            return candidate.max(1);
        }
    }
    magnitude as u64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_position() {
        assert_eq!(format_position(1_234_567), "1,234,567");
        assert_eq!(format_position(1000), "1,000");
        assert_eq!(format_position(500), "500");
    }

    #[test]
    fn test_nice_interval() {
        assert!(nice_interval(8.0) >= 8);
        assert!(nice_interval(100.0) >= 100);
    }
}
