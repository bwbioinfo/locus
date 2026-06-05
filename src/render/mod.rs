pub mod coverage;
pub mod features;
pub mod reads;
pub mod reference;
pub mod ruler;

/// Maps between genomic coordinates and terminal columns.
#[derive(Clone, Copy)]
pub struct ViewTransform {
    pub region_start: u64,
    pub region_end: u64,
    pub cols: u16,
}

impl ViewTransform {
    pub fn new(start: u64, end: u64, cols: u16) -> Self {
        Self {
            region_start: start,
            region_end: end,
            cols,
        }
    }

    pub fn bp_per_col(&self) -> f64 {
        if self.cols == 0 {
            return 1.0;
        }
        (self.region_end - self.region_start) as f64 / self.cols as f64
    }

    /// Convert a 0-based genomic position to a terminal column (0-indexed).
    /// Returns None if outside the visible range.
    #[allow(dead_code)]
    pub fn bp_to_col(&self, pos: u64) -> Option<u16> {
        if pos < self.region_start || pos >= self.region_end {
            return None;
        }
        let span = (self.region_end - self.region_start) as f64;
        if span == 0.0 {
            return None;
        }
        let frac = (pos - self.region_start) as f64 / span;
        Some((frac * self.cols as f64) as u16)
    }

    /// Convert a bp range to a column range (clamped to visible area).
    pub fn bp_range_to_cols(&self, start: u64, end: u64) -> (u16, u16) {
        let span = (self.region_end - self.region_start) as f64;
        if span == 0.0 {
            return (0, 0);
        }
        let col_start = if start <= self.region_start {
            0u16
        } else if start >= self.region_end {
            self.cols
        } else {
            ((start - self.region_start) as f64 / span * self.cols as f64) as u16
        };
        let col_end = if end >= self.region_end {
            self.cols
        } else if end <= self.region_start {
            0u16
        } else {
            ((end - self.region_start) as f64 / span * self.cols as f64) as u16
        };
        (
            col_start,
            col_end.max(col_start.saturating_add(1)).min(self.cols),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bp_to_col() {
        let t = ViewTransform::new(1000, 2000, 100);
        assert_eq!(t.bp_to_col(1000), Some(0));
        assert_eq!(t.bp_to_col(1500), Some(50));
        assert_eq!(t.bp_to_col(999), None);
        assert_eq!(t.bp_to_col(2000), None);
    }

    #[test]
    fn test_bp_range_to_cols() {
        let t = ViewTransform::new(0, 100, 10);
        let (s, e) = t.bp_range_to_cols(0, 50);
        assert_eq!(s, 0);
        assert_eq!(e, 5);
    }
}
