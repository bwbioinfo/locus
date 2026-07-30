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
    pub insertion_gap: Option<InsertionGap>,
    pub selection_bracket: Option<u64>,
    insertion_bracket_count: u8,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct InsertionGap {
    pub ref_pos: u64,
    pub len: u64,
}

impl InsertionGap {
    /// 0-based reference base preceding this insertion.
    ///
    /// `ref_pos` remains the first reference base after the insertion so the
    /// visual gap is placed between aligned bases.
    pub fn anchor_ref_pos(self) -> u64 {
        self.ref_pos.saturating_sub(1)
    }
}

impl ViewTransform {
    pub fn new(start: u64, end: u64, cols: u16) -> Self {
        Self {
            region_start: start,
            region_end: end,
            cols,
            insertion_gap: None,
            selection_bracket: None,
            insertion_bracket_count: 1,
        }
    }

    pub fn with_insertion_gap(mut self, insertion_gap: Option<InsertionGap>) -> Self {
        self.insertion_gap = insertion_gap;
        self
    }

    /// Reserve bracket columns around a selected reference base.
    pub fn with_selection_bracket(mut self, selection_bracket: Option<u64>) -> Self {
        self.selection_bracket = selection_bracket;
        self
    }

    /// Use an outer bracket pair for a selected expanded insertion.
    pub fn with_double_insertion_brackets(mut self, double_bracketed: bool) -> Self {
        self.insertion_bracket_count = if double_bracketed { 2 } else { 1 };
        self
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
        let col = self.base_col_for_ref(pos)?;
        Some(self.apply_visual_offsets(pos, col))
    }

    /// Convert a terminal column to its 0-based genomic position.
    ///
    /// Columns in an expanded insertion gap resolve to the insertion anchor.
    pub fn col_to_bp(&self, col: u16) -> Option<u64> {
        if col >= self.cols || self.region_start >= self.region_end || self.cols == 0 {
            return None;
        }

        let mut base_col = col;
        if let Some(gap) = self
            .insertion_gap
            .filter(|gap| (self.region_start..self.region_end).contains(&gap.ref_pos))
        {
            let gap_start = self.insertion_start_col(gap)?;
            let gap_len = self.insertion_visual_len(gap) as u16;
            let gap_end = gap_start.saturating_add(gap_len);
            if (gap_start..gap_end).contains(&col) {
                return Some(gap.anchor_ref_pos());
            }
            if col >= gap_end {
                base_col = col.saturating_sub(gap_len);
            }
        }

        if let Some(selection_bracket) = self.selection_bracket {
            let selected_col = self.selection_col_before_insertion(selection_bracket)?;
            let left_col = selected_col.saturating_sub(1);
            let right_col = selected_col.saturating_add(1);
            if base_col == left_col || base_col == right_col {
                return Some(selection_bracket);
            }
            if base_col > selected_col {
                base_col = base_col.saturating_sub(2);
            } else if base_col == selected_col {
                base_col = base_col.saturating_sub(1);
            }
        }

        let span = self.region_end - self.region_start;
        let offset = (u64::from(base_col) * span).div_ceil(u64::from(self.cols));
        Some((self.region_start + offset).min(self.region_end.saturating_sub(1)))
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
            let col = ((start - self.region_start) as f64 / span * self.cols as f64) as u16;
            self.apply_visual_offsets(start, col)
        };
        let col_end = if end >= self.region_end {
            self.cols
        } else if end <= self.region_start {
            0u16
        } else {
            let col = ((end - self.region_start) as f64 / span * self.cols as f64) as u16;
            self.apply_visual_offsets(end, col)
        };
        (
            col_start,
            col_end.max(col_start.saturating_add(1)).min(self.cols),
        )
    }

    pub fn insertion_col(&self, insertion_ref_pos: u64, insertion_offset: u64) -> Option<u16> {
        let gap = self.insertion_gap?;
        if gap.ref_pos != insertion_ref_pos || insertion_offset >= gap.len {
            return None;
        }
        let left_border_col = self.insertion_border_cols(insertion_ref_pos)?.0;
        let col = left_border_col
            .saturating_add(self.insertion_bracket_count())
            .saturating_add(insertion_offset as u16);
        (col < self.cols).then_some(col)
    }

    pub fn insertion_border_cols(&self, insertion_ref_pos: u64) -> Option<(u16, u16)> {
        let gap = self.insertion_gap?;
        if gap.ref_pos != insertion_ref_pos {
            return None;
        }
        let base_col = self.insertion_start_col(gap)?;
        let right_col = base_col
            .saturating_add(gap.len as u16)
            .saturating_add(self.insertion_bracket_count().saturating_mul(2))
            .saturating_sub(1);
        (right_col < self.cols).then_some((base_col, right_col))
    }

    /// Terminal columns containing the fluorescent brackets around a selected base.
    pub fn selection_bracket_cols(&self) -> Option<(u16, u16)> {
        let selection_bracket = self.selection_bracket?;
        let selected_col = self.bp_to_col(selection_bracket)?;
        let left_col = selected_col.checked_sub(1)?;
        let right_col = selected_col.checked_add(1)?;
        (right_col < self.cols).then_some((left_col, right_col))
    }

    pub fn insertion_bracket_count(&self) -> u16 {
        u16::from(self.insertion_bracket_count)
    }

    fn base_col_for_ref(&self, pos: u64) -> Option<u16> {
        if pos < self.region_start || pos >= self.region_end || self.cols == 0 {
            return None;
        }
        let span = (self.region_end - self.region_start) as f64;
        if span == 0.0 {
            return None;
        }
        let frac = (pos - self.region_start) as f64 / span;
        Some((frac * self.cols as f64) as u16)
    }

    fn apply_insertion_gap(&self, pos: u64, col: u16) -> u16 {
        match self.insertion_gap {
            Some(gap) if pos >= gap.ref_pos => col
                .saturating_add(self.insertion_visual_len(gap) as u16)
                .min(self.cols),
            _ => col,
        }
    }

    fn insertion_visual_len(&self, gap: InsertionGap) -> u64 {
        gap.len
            .saturating_add(u64::from(self.insertion_bracket_count).saturating_mul(2))
    }

    fn insertion_start_col(&self, gap: InsertionGap) -> Option<u16> {
        let col = self.base_col_for_ref(gap.ref_pos)?;
        Some(self.apply_selection_bracket(gap.ref_pos, col))
    }

    fn selection_col_before_insertion(&self, selection_bracket: u64) -> Option<u16> {
        let col = self.base_col_for_ref(selection_bracket)?;
        Some(self.apply_selection_bracket(selection_bracket, col))
    }

    fn apply_visual_offsets(&self, pos: u64, col: u16) -> u16 {
        self.apply_insertion_gap(pos, self.apply_selection_bracket(pos, col))
    }

    fn apply_selection_bracket(&self, pos: u64, col: u16) -> u16 {
        match self.selection_bracket {
            Some(selection_bracket) if pos > selection_bracket => col.saturating_add(2),
            Some(selection_bracket) if pos == selection_bracket => col.saturating_add(1),
            _ => col,
        }
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
    fn test_col_to_bp() {
        let t = ViewTransform::new(1000, 2000, 100);
        assert_eq!(t.col_to_bp(0), Some(1000));
        assert_eq!(t.col_to_bp(50), Some(1500));
        assert_eq!(t.col_to_bp(99), Some(1990));
        assert_eq!(t.col_to_bp(100), None);
    }

    #[test]
    fn col_to_bp_round_trips_to_the_rendered_column() {
        let t = ViewTransform::new(0, 100, 30);

        for col in 0..t.cols {
            let pos = t.col_to_bp(col).expect("visible position");
            assert_eq!(t.bp_to_col(pos), Some(col));
        }
    }

    #[test]
    fn insertion_gap_columns_resolve_to_the_anchor_position() {
        let t = ViewTransform::new(100, 200, 20).with_insertion_gap(Some(InsertionGap {
            ref_pos: 150,
            len: 2,
        }));

        assert_eq!(t.col_to_bp(10), Some(149));
        assert_eq!(t.col_to_bp(11), Some(149));
        assert_eq!(t.col_to_bp(12), Some(149));
        assert_eq!(t.col_to_bp(13), Some(149));
        assert_eq!(t.col_to_bp(14), Some(150));
        assert_eq!(t.col_to_bp(15), Some(155));
    }

    #[test]
    fn selection_brackets_reserve_columns_without_replacing_the_selected_base() {
        let t = ViewTransform::new(100, 110, 10).with_selection_bracket(Some(104));

        assert_eq!(t.bp_to_col(103), Some(3));
        assert_eq!(t.bp_to_col(104), Some(5));
        assert_eq!(t.bp_to_col(105), Some(7));
        assert_eq!(t.selection_bracket_cols(), Some((4, 6)));
        assert_eq!(t.col_to_bp(4), Some(104));
        assert_eq!(t.col_to_bp(5), Some(104));
        assert_eq!(t.col_to_bp(6), Some(104));
        assert_eq!(t.col_to_bp(7), Some(105));
    }

    #[test]
    fn selected_insertions_reserve_an_outer_bracket_pair() {
        let t = ViewTransform::new(100, 120, 20)
            .with_insertion_gap(Some(InsertionGap {
                ref_pos: 110,
                len: 2,
            }))
            .with_double_insertion_brackets(true);

        assert_eq!(t.insertion_border_cols(110), Some((10, 15)));
        assert_eq!(t.insertion_col(110, 0), Some(12));
        assert_eq!(t.insertion_col(110, 1), Some(13));
        assert_eq!(t.col_to_bp(10), Some(109));
        assert_eq!(t.col_to_bp(15), Some(109));
        assert_eq!(t.bp_to_col(110), Some(16));
    }

    #[test]
    fn selection_brackets_and_insertion_gaps_keep_regular_columns_clickable() {
        let t = ViewTransform::new(100, 120, 20)
            .with_insertion_gap(Some(InsertionGap {
                ref_pos: 110,
                len: 1,
            }))
            .with_selection_bracket(Some(105));

        assert_eq!(t.bp_to_col(105), Some(6));
        assert_eq!(t.bp_to_col(106), Some(8));
        assert_eq!(t.bp_to_col(110), Some(15));
        assert_eq!(t.col_to_bp(5), Some(105));
        assert_eq!(t.col_to_bp(7), Some(105));
        assert_eq!(t.col_to_bp(8), Some(106));
        assert_eq!(t.col_to_bp(12), Some(109));
        assert_eq!(t.col_to_bp(15), Some(110));
    }

    #[test]
    fn test_bp_range_to_cols() {
        let t = ViewTransform::new(0, 100, 10);
        let (s, e) = t.bp_range_to_cols(0, 50);
        assert_eq!(s, 0);
        assert_eq!(e, 5);
    }
}
