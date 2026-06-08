use anyhow::Result;

use crate::{
    bam::BamSource,
    cache::RegionCache,
    gff::GffStore,
    reference::ReferenceStore,
    region::{Region, parse_region},
    render::{
        InsertionGap, ViewTransform,
        reads::{selected_insertion_gap, visible_insertion_gaps},
    },
    screenshot,
};

/// UI interaction mode.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Mode {
    Normal,
    /// Typing a "go to region" string
    GoTo,
    /// Typing a gene/feature name to search
    FeatureSearch,
    /// Choosing a contig from a list
    ContigSelect,
    Help,
}

/// Zoom / bp-per-terminal-column.
const MIN_BP_PER_COL: f64 = 1.0;
const MAX_BP_PER_COL: f64 = 100_000.0;

/// Application state.
pub struct App {
    pub source: BamSource,
    pub cache: RegionCache,
    pub gff: Option<GffStore>,
    pub reference: Option<ReferenceStore>,

    // --- browser state ---
    pub contig_idx: usize,
    /// 0-based start of the visible window
    pub view_start: u64,
    /// 0-based exclusive end of the visible window
    pub view_end: u64,
    /// bp represented by each terminal column (zoom)
    pub bp_per_col: f64,

    // --- ui state ---
    pub mode: Mode,
    pub command_buffer: String,
    pub show_help: bool,
    pub status_msg: Option<String>,
    pub terminal_cols: u16,
    pub terminal_rows: u16,
    pub expand_insertions: bool,
    pub selected_insertion_ref_pos: Option<u64>,

    /// Set to true to request a clean exit
    pub should_quit: bool,
    /// True when the cache needs to be refreshed before rendering
    pub needs_fetch: bool,

    // --- feature search state ---
    /// Indices into gff.features matching the current search query
    pub feature_matches: Vec<usize>,
    /// Which match is currently highlighted / will be jumped to
    pub feature_match_idx: usize,
}

impl App {
    pub fn new(
        source: BamSource,
        gff: Option<GffStore>,
        reference: Option<ReferenceStore>,
        initial_region: Option<Region>,
    ) -> Result<Self> {
        let (view_start, view_end) = if let Some(ref r) = initial_region {
            (r.start, r.end)
        } else if let Some(contig) = source.contigs.first() {
            let len = contig.length;
            (0, 1_000.min(len))
        } else {
            (0, 1_000)
        };

        let contig_idx = if let Some(ref r) = initial_region {
            source
                .contigs
                .iter()
                .position(|c| c.name == r.contig)
                .unwrap_or(0)
        } else {
            0
        };

        let span = (view_end - view_start) as f64;
        let bp_per_col = (span / 80.0).max(MIN_BP_PER_COL);

        let mut app = Self {
            source,
            cache: RegionCache::default(),
            gff,
            reference,
            contig_idx,
            view_start,
            view_end,
            bp_per_col,
            mode: Mode::Normal,
            command_buffer: String::new(),
            show_help: false,
            status_msg: None,
            terminal_cols: 80,
            terminal_rows: 24,
            expand_insertions: false,
            selected_insertion_ref_pos: None,
            should_quit: false,
            needs_fetch: true,
            feature_matches: Vec::new(),
            feature_match_idx: 0,
        };

        app.clamp_view();
        Ok(app)
    }

    // ─── Accessors ────────────────────────────────────────────────────────────

    pub fn current_contig(&self) -> &str {
        self.source
            .contigs
            .get(self.contig_idx)
            .map(|c| c.name.as_str())
            .unwrap_or("?")
    }

    pub fn current_contig_len(&self) -> u64 {
        self.source
            .contigs
            .get(self.contig_idx)
            .map(|c| c.length)
            .unwrap_or(0)
    }

    pub fn current_region(&self) -> Region {
        Region::new(self.current_contig(), self.view_start, self.view_end)
    }

    pub fn view_span(&self) -> u64 {
        self.view_end.saturating_sub(self.view_start)
    }

    pub fn view_cols(&self) -> usize {
        self.terminal_cols.saturating_sub(2) as usize
    }

    pub fn toggle_insertions(&mut self) {
        self.expand_insertions = !self.expand_insertions;
        if !self.expand_insertions {
            self.selected_insertion_ref_pos = None;
        }
        self.status_msg = Some(if self.expand_insertions {
            "insertions expanded".to_string()
        } else {
            "insertions collapsed".to_string()
        });
    }

    pub fn cycle_insertion_expansion(&mut self, forward: bool) {
        if !self.expand_insertions {
            self.expand_insertions = true;
        }

        let transform = self.base_view_transform();
        let gaps = visible_insertion_gaps(&self.cache.reads, &self.cache.pileup_rows, &transform);
        if gaps.is_empty() {
            self.selected_insertion_ref_pos = None;
            self.status_msg = Some("no visible insertions to expand".to_string());
            return;
        }

        let current_idx = self
            .selected_insertion_ref_pos
            .and_then(|pos| gaps.iter().position(|gap| gap.ref_pos == pos));
        let next_idx = match (current_idx, forward) {
            (Some(idx), true) => (idx + 1) % gaps.len(),
            (Some(0), false) => gaps.len() - 1,
            (Some(idx), false) => idx - 1,
            (None, true) => 0,
            (None, false) => gaps.len() - 1,
        };
        let gap = gaps[next_idx];
        self.selected_insertion_ref_pos = Some(gap.ref_pos);
        self.status_msg = Some(format!(
            "expanded insertion {} bp at {}:{}",
            gap.len,
            self.current_contig(),
            gap.ref_pos + 1
        ));
    }

    pub fn selected_insertion_gap(&self, transform: &ViewTransform) -> Option<InsertionGap> {
        if !self.expand_insertions {
            return None;
        }
        let gaps = visible_insertion_gaps(&self.cache.reads, &self.cache.pileup_rows, transform);
        if let Some(selected_ref_pos) = self.selected_insertion_ref_pos {
            if let Some(gap) = gaps
                .iter()
                .copied()
                .find(|gap| gap.ref_pos == selected_ref_pos)
            {
                return Some(gap);
            }
        }
        selected_insertion_gap(&self.cache.reads, &self.cache.pileup_rows, transform)
    }

    fn base_view_transform(&self) -> ViewTransform {
        ViewTransform::new(
            self.view_start,
            self.view_end,
            self.terminal_cols.saturating_sub(2),
        )
    }

    // ─── Navigation ───────────────────────────────────────────────────────────

    pub fn pan(&mut self, delta_bp: i64) {
        let len = self.current_contig_len();
        let span = self.view_span();
        if delta_bp > 0 {
            let d = delta_bp as u64;
            self.view_start = (self.view_start + d).min(len.saturating_sub(span));
        } else {
            let d = (-delta_bp) as u64;
            self.view_start = self.view_start.saturating_sub(d);
        }
        self.view_end = (self.view_start + span).min(len);
        self.mark_dirty();
    }

    pub fn zoom_in(&mut self) {
        let center = (self.view_start + self.view_end) / 2;
        let half_span = (self.view_span() / 4).max(50);
        self.view_start = center.saturating_sub(half_span);
        self.view_end = center + half_span;
        self.bp_per_col = (self.view_span() as f64 / self.view_cols() as f64).max(MIN_BP_PER_COL);
        self.clamp_view();
        self.mark_dirty();
    }

    pub fn zoom_out(&mut self) {
        let center = (self.view_start + self.view_end) / 2;
        let half_span = (self.view_span()).min(MAX_BP_PER_COL as u64 * self.view_cols() as u64 / 2);
        let new_half = (half_span * 2).min(self.current_contig_len());
        self.view_start = center.saturating_sub(new_half / 2);
        self.view_end = center + new_half / 2;
        self.bp_per_col = (self.view_span() as f64 / self.view_cols() as f64).max(MIN_BP_PER_COL);
        self.clamp_view();
        self.mark_dirty();
    }

    /// If the new view is within the cached padded region, just re-layout without disk IO.
    /// Only set needs_fetch=true when the view has drifted outside the loaded window.
    fn mark_dirty(&mut self) {
        let within_cache = self.cache.loaded_region.as_ref().map_or(false, |loaded| {
            loaded.contig == self.current_contig()
                && self.view_start >= loaded.start
                && self.view_end <= loaded.end
        });
        if within_cache {
            self.relayout();
        } else {
            self.needs_fetch = true;
        }
    }

    /// Re-layout pileup and coverage from the already-loaded reads (no disk IO).
    pub fn relayout(&mut self) {
        let visible = self.current_region();
        let reference_rows = usize::from(self.reference.is_some());
        let max_rows = self
            .terminal_rows
            .saturating_sub(12 + reference_rows as u16) as usize;
        let cols = self.view_cols();
        self.cache.layout_pileup(&visible, max_rows.max(1));
        self.cache.compute_coverage(&visible, cols.max(1));
    }

    pub fn jump_to_region(&mut self, region: &Region) -> Result<()> {
        let idx = self
            .source
            .contigs
            .iter()
            .position(|c| c.name == region.contig)
            .ok_or_else(|| crate::error::LocusError::UnknownContig(region.contig.clone()))?;

        self.contig_idx = idx;
        let len = self.current_contig_len();

        self.view_start = region.start.min(len.saturating_sub(1));
        self.view_end = if region.end == u64::MAX {
            (self.view_start + 1_000).min(len)
        } else {
            region.end.min(len)
        };

        self.bp_per_col = (self.view_span() as f64 / self.view_cols() as f64).max(MIN_BP_PER_COL);
        self.needs_fetch = true;
        Ok(())
    }

    pub fn select_contig(&mut self, idx: usize) {
        if idx < self.source.contigs.len() {
            self.contig_idx = idx;
            let len = self.current_contig_len();
            self.view_start = 0;
            self.view_end = 1_000.min(len);
            self.bp_per_col =
                (self.view_span() as f64 / self.view_cols() as f64).max(MIN_BP_PER_COL);
            self.needs_fetch = true;
        }
    }

    // ─── Feature search ───────────────────────────────────────────────────────

    /// Run a search against the GFF store, updating `feature_matches`.
    pub fn run_feature_search(&mut self) {
        let query = self.command_buffer.trim().to_string();
        if let Some(ref gff) = self.gff {
            self.feature_matches = gff.search(&query);
        } else {
            self.feature_matches.clear();
        }
        self.feature_match_idx = 0;
    }

    /// Jump to the currently selected feature match.
    pub fn jump_to_current_match(&mut self) -> Result<()> {
        let idx = match self.feature_matches.get(self.feature_match_idx) {
            Some(&i) => i,
            None => {
                self.status_msg = Some("No matching features".to_string());
                return Ok(());
            }
        };

        let region = {
            let gff = self.gff.as_ref().unwrap();
            let feat = &gff.features[idx];
            // Add 10 % padding each side
            let pad = (feat.end - feat.start) / 10 + 1;
            let padded_start = feat.start.saturating_sub(pad);
            let padded_end = feat.end + pad;
            Region::new(feat.seqname.clone(), padded_start, padded_end)
        };

        self.jump_to_region(&region).map_err(|e| {
            self.status_msg = Some(format!("{e}"));
            e
        })?;
        Ok(())
    }

    /// Cycle to the next search result and jump to it.
    pub fn next_feature_match(&mut self) -> Result<()> {
        if !self.feature_matches.is_empty() {
            self.feature_match_idx = (self.feature_match_idx + 1) % self.feature_matches.len();
            self.jump_to_current_match()?;
        }
        Ok(())
    }

    /// Cycle to the previous search result and jump to it.
    pub fn prev_feature_match(&mut self) -> Result<()> {
        if !self.feature_matches.is_empty() {
            self.feature_match_idx = self
                .feature_match_idx
                .checked_sub(1)
                .unwrap_or(self.feature_matches.len() - 1);
            self.jump_to_current_match()?;
        }
        Ok(())
    }

    // ─── Data fetching ────────────────────────────────────────────────────────

    pub fn refresh(&mut self) -> Result<()> {
        let visible = self.current_region();
        let len = self.current_contig_len();
        let pad = self.view_span() / 2;
        let padded = visible.padded(pad, len);

        let reads = self.source.fetch_reads(&padded).map_err(|e| {
            self.status_msg = Some(format!("Error: {e}"));
            e
        })?;

        let reference_rows = usize::from(self.reference.is_some());
        let max_pileup_rows = self
            .terminal_rows
            .saturating_sub(12 + reference_rows as u16) as usize;
        let view_cols = self.view_cols();

        self.cache.reads = reads;
        self.cache.reference = if let Some(reference) = self.reference.as_ref() {
            reference.fetch(&padded)?
        } else {
            None
        };
        self.cache.loaded_region = Some(padded);
        self.cache.layout_pileup(&visible, max_pileup_rows.max(1));
        self.cache.compute_coverage(&visible, view_cols.max(1));

        self.needs_fetch = false;
        self.status_msg = None;
        Ok(())
    }

    pub fn save_screenshot(&mut self) {
        match screenshot::save(self) {
            Ok(paths) => {
                self.status_msg = Some(format!(
                    "screenshot: {} + {}",
                    paths.text.display(),
                    paths.html.display()
                ));
            }
            Err(e) => {
                self.status_msg = Some(format!("screenshot failed: {e}"));
            }
        }
    }

    // ─── Input handling ───────────────────────────────────────────────────────

    pub fn handle_goto_input(&mut self, c: char) {
        self.command_buffer.push(c);
    }

    pub fn confirm_goto(&mut self) -> Result<()> {
        let input = self.command_buffer.trim().to_string();
        self.command_buffer.clear();
        self.mode = Mode::Normal;

        let region = parse_region(&input).map_err(|e| {
            self.status_msg = Some(format!("{e}"));
            e
        })?;
        self.jump_to_region(&region).map_err(|e| {
            self.status_msg = Some(format!("{e}"));
            e
        })?;
        Ok(())
    }

    pub fn cancel_input(&mut self) {
        self.command_buffer.clear();
        self.mode = Mode::Normal;
    }

    // ─── Helpers ──────────────────────────────────────────────────────────────

    fn clamp_view(&mut self) {
        let len = self.current_contig_len();
        if len == 0 {
            return;
        }
        let span = self.view_span().max(1);
        self.view_start = self.view_start.min(len.saturating_sub(1));
        self.view_end = (self.view_start + span).min(len);
        if self.view_end == self.view_start {
            self.view_end = self.view_start + 1;
        }
    }
}
