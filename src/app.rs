use anyhow::Result;

use crate::{
    bam::BamSource,
    cache::{PhasePositionAlleleTallies, PositionAlleleTally, RegionCache},
    gff::GffStore,
    reference::ReferenceStore,
    region::{Region, parse_region},
    render::{
        InsertionGap, ViewTransform,
        reads::{selected_insertion_gap, visible_insertion_gaps},
    },
    screenshot,
    theme::Theme,
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
    /// Typing a minimum mapping-quality threshold
    MapqFilter,
    Help,
}

/// A vertically scrollable read section.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum ReadTrack {
    #[default]
    Combined,
    Hp1,
    Hp2,
    Unphased,
}

impl ReadTrack {
    const ALL: [Self; 4] = [Self::Combined, Self::Hp1, Self::Hp2, Self::Unphased];

    fn index(self) -> usize {
        match self {
            Self::Combined => 0,
            Self::Hp1 => 1,
            Self::Hp2 => 2,
            Self::Unphased => 3,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::Combined => "combined",
            Self::Hp1 => "HP1",
            Self::Hp2 => "HP2",
            Self::Unphased => "unphased",
        }
    }
}

/// Zoom / bp-per-terminal-column.
const MIN_BP_PER_COL: f64 = 1.0;
const MAX_BP_PER_COL: f64 = 100_000.0;
pub const LARGE_PAN_BP: i64 = 1_000;

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
    /// 0-based base preceding the insertion selected by Tab navigation.
    pub selected_insertion_ref_pos: Option<u64>,
    /// 0-based genomic position selected with the mouse.
    pub selected_ref_pos: Option<u64>,
    /// Alleles observed at the selected reference position under the active MAPQ filter.
    pub selected_allele_tally: Option<PositionAlleleTally>,
    /// Phase-separated alleles observed at the selected reference position.
    pub selected_phase_allele_tallies: Option<PhasePositionAlleleTallies>,
    pub show_selection_brackets: bool,
    pub show_methylation: bool,
    pub show_phasing: bool,
    pub active_read_track: ReadTrack,
    read_scroll_offsets: [usize; 4],
    pub theme: Theme,
    pub min_mapq: u8,

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
        theme: Theme,
        min_mapq: u8,
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
            selected_ref_pos: None,
            selected_allele_tally: None,
            selected_phase_allele_tallies: None,
            show_selection_brackets: true,
            show_methylation: false,
            show_phasing: false,
            active_read_track: ReadTrack::Combined,
            read_scroll_offsets: [0; 4],
            theme,
            min_mapq,
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
        } else if let Some(position) = self.selected_ref_pos {
            self.activate_insertion_at_selected_position(position);
        }
        self.status_msg = Some(if self.expand_insertions {
            "insertions expanded".to_string()
        } else {
            "insertions collapsed".to_string()
        });
    }

    pub fn toggle_methylation(&mut self) {
        self.show_methylation = !self.show_methylation;
        self.status_msg = Some(if self.show_methylation {
            "methylation shown".to_string()
        } else {
            "methylation hidden".to_string()
        });
    }

    pub fn toggle_phasing(&mut self) {
        self.show_phasing = !self.show_phasing;
        self.active_read_track = if self.show_phasing {
            ReadTrack::Hp1
        } else {
            ReadTrack::Combined
        };
        self.relayout();
        self.status_msg = Some(if self.show_phasing {
            "phasing split into HP1 and HP2 tracks; unphased reads separated".to_string()
        } else {
            "phasing hidden".to_string()
        });
    }

    pub fn toggle_theme(&mut self) {
        self.theme = self.theme.toggle();
        self.status_msg = Some(format!("{} theme", self.theme.name()));
    }

    pub fn toggle_selection_brackets(&mut self) {
        self.show_selection_brackets = !self.show_selection_brackets;
        self.status_msg = Some(if self.show_selection_brackets {
            "selection brackets shown".to_string()
        } else {
            "selection brackets hidden".to_string()
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
            .and_then(|pos| gaps.iter().position(|gap| gap.anchor_ref_pos() == pos));
        let next_idx = match (current_idx, forward) {
            (Some(idx), true) => (idx + 1) % gaps.len(),
            (Some(0), false) => gaps.len() - 1,
            (Some(idx), false) => idx - 1,
            (None, true) => 0,
            (None, false) => gaps.len() - 1,
        };
        let gap = gaps[next_idx];
        self.selected_insertion_ref_pos = Some(gap.anchor_ref_pos());
        self.status_msg = Some(format!(
            "expanded insertion {} bp at {}:{}",
            gap.len,
            self.current_contig(),
            gap.anchor_ref_pos() + 1
        ));
    }

    pub fn selected_insertion_gap(&self, transform: &ViewTransform) -> Option<InsertionGap> {
        if !self.expand_insertions {
            return None;
        }
        let gaps = visible_insertion_gaps(&self.cache.reads, &self.cache.pileup_rows, transform);
        if let Some(selected_ref_pos) = self.selected_insertion_ref_pos
            && let Some(gap) = gaps
                .iter()
                .copied()
                .find(|gap| gap.anchor_ref_pos() == selected_ref_pos)
        {
            return Some(gap);
        }
        selected_insertion_gap(&self.cache.reads, &self.cache.pileup_rows, transform)
    }

    pub fn select_reference_position(&mut self, pos: u64) {
        if (self.view_start..self.view_end).contains(&pos) {
            self.selected_ref_pos = Some(pos);
            self.activate_insertion_at_selected_position(pos);
            self.refresh_selected_allele_tally();
            self.status_msg = None;
        }
    }

    /// Clear the selected reference position, its tallies, and any selected insertion anchor.
    pub fn clear_selected_position(&mut self) {
        self.clear_selected_reference_position();
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
        if self.selected_ref_pos.is_some() {
            self.move_selected_reference_position(delta_bp.signum());
            return;
        }
        self.pan_view(delta_bp);
    }

    pub fn pan_large(&mut self, direction: i64) {
        self.pan_view(direction.signum() * LARGE_PAN_BP);
    }

    fn pan_view(&mut self, delta_bp: i64) {
        self.shift_view(delta_bp);
        self.mark_dirty(true);
    }

    fn shift_view(&mut self, delta_bp: i64) {
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
    }

    fn move_selected_reference_position(&mut self, delta_bp: i64) {
        let Some(selected) = self.selected_ref_pos else {
            return;
        };
        let next = if delta_bp > 0 {
            selected
                .saturating_add(1)
                .min(self.current_contig_len().saturating_sub(1))
        } else if delta_bp < 0 {
            selected.saturating_sub(1)
        } else {
            selected
        };
        if next == selected {
            return;
        }

        if next < self.view_start {
            self.shift_view(-1);
        } else if next >= self.view_end {
            self.shift_view(1);
        }
        self.selected_ref_pos = Some(next);
        self.mark_dirty(false);
        self.activate_insertion_at_selected_position(next);
    }

    /// Make an insertion at the selected anchor the active shared expansion.
    fn activate_insertion_at_selected_position(&mut self, position: u64) {
        if !self.expand_insertions {
            return;
        }

        let transform = self.base_view_transform();
        if visible_insertion_gaps(&self.cache.reads, &self.cache.pileup_rows, &transform)
            .iter()
            .any(|gap| gap.anchor_ref_pos() == position)
        {
            self.selected_insertion_ref_pos = Some(position);
        }
    }

    pub fn zoom_in(&mut self) {
        let center = self.zoom_center();
        let half_span = (self.view_span() / 4).max(50);
        self.view_start = center.saturating_sub(half_span);
        self.view_end = center + half_span;
        self.bp_per_col = (self.view_span() as f64 / self.view_cols() as f64).max(MIN_BP_PER_COL);
        self.clamp_view();
        self.mark_dirty(false);
    }

    pub fn zoom_out(&mut self) {
        let center = self.zoom_center();
        let half_span = (self.view_span()).min(MAX_BP_PER_COL as u64 * self.view_cols() as u64 / 2);
        let new_half = (half_span * 2).min(self.current_contig_len());
        self.view_start = center.saturating_sub(new_half / 2);
        self.view_end = center + new_half / 2;
        self.bp_per_col = (self.view_span() as f64 / self.view_cols() as f64).max(MIN_BP_PER_COL);
        self.clamp_view();
        self.mark_dirty(false);
    }

    fn zoom_center(&self) -> u64 {
        self.selected_ref_pos
            .unwrap_or_else(|| (self.view_start + self.view_end) / 2)
    }

    /// If the new view is within the cached padded region, just re-layout without disk IO.
    /// Only set needs_fetch=true when the view has drifted outside the loaded window.
    fn mark_dirty(&mut self, clear_selection: bool) {
        if clear_selection {
            self.clear_selected_reference_position();
        }
        let within_cache = self.cache.loaded_region.as_ref().is_some_and(|loaded| {
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
        let max_rows = crate::ui::available_read_rows(
            self.terminal_rows,
            self.reference.is_some(),
            self.gff.is_some(),
        );
        let cols = self.view_cols();
        self.cache
            .layout_pileup(&visible, max_rows.max(1), self.min_mapq, self.show_phasing);
        self.cache
            .compute_coverage(&visible, cols.max(1), self.min_mapq);
        self.clamp_read_scrolls();
        self.refresh_selected_allele_tally();
    }

    pub fn read_track_scroll(&self, track: ReadTrack) -> usize {
        self.read_scroll_offsets[track.index()]
    }

    pub fn scroll_read_track(&mut self, track: ReadTrack, delta_rows: i32) {
        let track = self.normalized_read_track(track);
        self.active_read_track = track;
        let row_count = self.read_track_row_count(track);
        let viewport_rows = self.read_track_viewport_rows(track).max(1);
        let max_scroll = row_count.saturating_sub(viewport_rows);
        let current_offset = self.read_scroll_offsets[track.index()];
        let offset = if delta_rows >= 0 {
            current_offset
                .saturating_add(delta_rows as usize)
                .min(max_scroll)
        } else {
            current_offset.saturating_sub(delta_rows.unsigned_abs() as usize)
        };
        self.read_scroll_offsets[track.index()] = offset;

        if row_count > 0 {
            let first = offset + 1;
            let last = (offset + viewport_rows).min(row_count);
            self.status_msg = Some(format!(
                "{} read rows {first}-{last} of {row_count}",
                track.name()
            ));
        }
    }

    pub fn scroll_active_read_track(&mut self, delta_rows: i32) {
        self.scroll_read_track(self.active_read_track, delta_rows);
    }

    pub fn cycle_active_read_track(&mut self, forward: bool) {
        let tracks = self.available_read_tracks();
        let current = tracks
            .iter()
            .position(|&track| track == self.active_read_track)
            .unwrap_or(0);
        let next = if forward {
            (current + 1) % tracks.len()
        } else {
            current.checked_sub(1).unwrap_or(tracks.len() - 1)
        };
        self.active_read_track = tracks[next];
        self.status_msg = Some(format!(
            "active read track: {}",
            self.active_read_track.name()
        ));
    }

    fn available_read_tracks(&self) -> Vec<ReadTrack> {
        let Some(layout) = self
            .cache
            .phase_layout
            .as_ref()
            .filter(|_| self.show_phasing)
        else {
            return vec![ReadTrack::Combined];
        };

        let mut tracks = vec![ReadTrack::Hp1, ReadTrack::Hp2];
        if layout.unphased_rows.is_some() {
            tracks.push(ReadTrack::Unphased);
        }
        tracks
    }

    fn normalized_read_track(&self, track: ReadTrack) -> ReadTrack {
        self.available_read_tracks()
            .into_iter()
            .find(|&available| available == track)
            .unwrap_or(ReadTrack::Combined)
    }

    fn read_track_row_count(&self, track: ReadTrack) -> usize {
        let track = self.normalized_read_track(track);
        match track {
            ReadTrack::Combined => self.cache.pileup_rows.len(),
            ReadTrack::Hp1 => self
                .cache
                .phase_layout
                .as_ref()
                .map_or(0, |layout| layout.hp1_rows.len()),
            ReadTrack::Hp2 => self
                .cache
                .phase_layout
                .as_ref()
                .map_or(0, |layout| layout.hp2_rows.len()),
            ReadTrack::Unphased => self
                .cache
                .phase_layout
                .as_ref()
                .and_then(|layout| layout.unphased_rows.as_ref())
                .map_or(0, std::ops::Range::len),
        }
    }

    fn read_track_viewport_rows(&self, track: ReadTrack) -> usize {
        let track = self.normalized_read_track(track);
        match track {
            ReadTrack::Combined => crate::ui::available_read_rows(
                self.terminal_rows,
                self.reference.is_some(),
                self.gff.is_some(),
            ),
            ReadTrack::Hp1 => self
                .cache
                .phase_layout
                .as_ref()
                .map_or(0, |layout| layout.hp1_viewport_rows),
            ReadTrack::Hp2 => self
                .cache
                .phase_layout
                .as_ref()
                .map_or(0, |layout| layout.hp2_viewport_rows),
            ReadTrack::Unphased => self
                .cache
                .phase_layout
                .as_ref()
                .map_or(0, |layout| layout.unphased_viewport_rows),
        }
    }

    fn clamp_read_scrolls(&mut self) {
        for track in ReadTrack::ALL {
            let max_scroll = self
                .read_track_row_count(track)
                .saturating_sub(self.read_track_viewport_rows(track).max(1));
            self.read_scroll_offsets[track.index()] =
                self.read_scroll_offsets[track.index()].min(max_scroll);
        }
        self.active_read_track = self.normalized_read_track(self.active_read_track);
    }

    pub fn jump_to_region(&mut self, region: &Region) -> Result<()> {
        let idx = self
            .source
            .contigs
            .iter()
            .position(|c| c.name == region.contig)
            .ok_or_else(|| crate::error::LocusError::UnknownContig(region.contig.clone()))?;

        self.contig_idx = idx;
        self.clear_selected_reference_position();
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
            self.clear_selected_reference_position();
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

        let region = self.gff.as_ref().unwrap().features[idx].padded_region();

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

        let max_pileup_rows = crate::ui::available_read_rows(
            self.terminal_rows,
            self.reference.is_some(),
            self.gff.is_some(),
        );
        let view_cols = self.view_cols();

        self.cache.reads = reads;
        self.cache.reference = if let Some(reference) = self.reference.as_ref() {
            reference.fetch(&padded)?
        } else {
            None
        };
        self.cache.loaded_region = Some(padded);
        self.cache.layout_pileup(
            &visible,
            max_pileup_rows.max(1),
            self.min_mapq,
            self.show_phasing,
        );
        self.cache
            .compute_coverage(&visible, view_cols.max(1), self.min_mapq);
        self.clamp_read_scrolls();
        self.refresh_selected_allele_tally();

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

    pub fn begin_mapq_filter(&mut self) {
        self.mode = Mode::MapqFilter;
        self.command_buffer.clear();
        self.status_msg = None;
    }

    pub fn confirm_mapq_filter(&mut self) {
        let Ok(min_mapq) = self.command_buffer.trim().parse::<u8>() else {
            self.status_msg = Some("minimum MAPQ must be between 0 and 255".to_string());
            return;
        };

        self.min_mapq = min_mapq;
        self.command_buffer.clear();
        self.mode = Mode::Normal;
        self.relayout();
        self.status_msg = Some(if min_mapq == 0 {
            "MAPQ filter disabled".to_string()
        } else {
            format!("minimum MAPQ set to {min_mapq}")
        });
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

    fn clear_selected_reference_position(&mut self) {
        self.selected_ref_pos = None;
        self.selected_insertion_ref_pos = None;
        self.selected_allele_tally = None;
        self.selected_phase_allele_tallies = None;
    }

    fn refresh_selected_allele_tally(&mut self) {
        let Some(position) = self.selected_ref_pos else {
            self.selected_allele_tally = None;
            self.selected_phase_allele_tallies = None;
            return;
        };

        self.selected_allele_tally = Some(self.cache.allele_tally_at(position, self.min_mapq));
        self.selected_phase_allele_tallies = self
            .show_phasing
            .then(|| self.cache.phase_allele_tallies_at(position, self.min_mapq));
    }
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use super::*;

    fn demo_app(min_mapq: u8) -> App {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/demo/demo.sorted.bam");
        let source = BamSource::open(path).expect("open demo BAM");

        App::new(source, None, None, None, Theme::Dark, min_mapq).expect("create app")
    }

    fn read_with_insertion(name: &str, matched_before_insertion: u64) -> crate::cache::RenderRead {
        crate::cache::RenderRead {
            name: name.to_string(),
            start: 0,
            end: 20,
            strand: crate::cache::Strand::Forward,
            mapq: 60,
            cigar_ops: vec![
                crate::cache::CigarOp::Match(matched_before_insertion),
                crate::cache::CigarOp::Insertion(2),
                crate::cache::CigarOp::Match(20 - matched_before_insertion),
            ],
            sequence: vec![b'A'; 22],
            methylation: Vec::new(),
            deleted_reference_sequences: Vec::new(),
            phase: Default::default(),
            is_secondary: false,
            is_supplementary: false,
            is_duplicate: false,
        }
    }

    #[test]
    fn mapq_prompt_applies_without_refetching() {
        let mut app = demo_app(0);
        app.needs_fetch = false;
        app.begin_mapq_filter();
        app.command_buffer.push_str("30");

        app.confirm_mapq_filter();

        assert_eq!(app.min_mapq, 30);
        assert_eq!(app.mode, Mode::Normal);
        assert!(!app.needs_fetch);
        assert_eq!(app.status_msg.as_deref(), Some("minimum MAPQ set to 30"));
    }

    #[test]
    fn mapq_prompt_rejects_out_of_range_value() {
        let mut app = demo_app(10);
        app.needs_fetch = false;
        app.begin_mapq_filter();
        app.command_buffer.push_str("256");

        app.confirm_mapq_filter();

        assert_eq!(app.min_mapq, 10);
        assert_eq!(app.mode, Mode::MapqFilter);
        assert!(!app.needs_fetch);
        assert_eq!(
            app.status_msg.as_deref(),
            Some("minimum MAPQ must be between 0 and 255")
        );
    }

    #[test]
    fn phasing_is_off_by_default_and_toggles_without_refetching() {
        let mut app = demo_app(0);
        assert!(!app.show_phasing);
        app.refresh().expect("load demo reads");
        assert!(app.cache.phase_layout.is_none());
        let combined_rows = app.cache.pileup_rows.clone();

        app.toggle_phasing();

        assert!(app.show_phasing);
        assert!(!app.needs_fetch);
        assert!(app.cache.phase_layout.is_some());
        assert_eq!(
            app.status_msg.as_deref(),
            Some("phasing split into HP1 and HP2 tracks; unphased reads separated")
        );

        app.toggle_phasing();

        assert!(!app.show_phasing);
        assert!(!app.needs_fetch);
        assert!(app.cache.phase_layout.is_none());
        assert_eq!(app.cache.pileup_rows, combined_rows);
        assert_eq!(app.status_msg.as_deref(), Some("phasing hidden"));
    }

    #[test]
    fn selected_tally_switches_to_phase_groups_with_phasing_enabled() {
        let mut app = demo_app(0);
        app.view_start = 0;
        app.view_end = 10;
        let mut hp1 = read_with_insertion("hp1", 1);
        hp1.phase.haplotype = Some(1);
        app.cache.reads = vec![hp1];

        app.select_reference_position(0);
        assert!(app.selected_allele_tally.is_some());
        assert!(app.selected_phase_allele_tallies.is_none());

        app.show_phasing = true;
        app.refresh_selected_allele_tally();

        assert_eq!(
            app.selected_phase_allele_tallies
                .as_ref()
                .and_then(|tallies| tallies.hp1.base_counts.get(&b'A')),
            Some(&1)
        );
    }

    #[test]
    fn selection_brackets_are_enabled_by_default_and_toggle_without_refetching() {
        let mut app = demo_app(0);
        app.needs_fetch = false;

        assert!(app.show_selection_brackets);
        app.toggle_selection_brackets();

        assert!(!app.show_selection_brackets);
        assert!(!app.needs_fetch);
        assert_eq!(app.status_msg.as_deref(), Some("selection brackets hidden"));
    }

    #[test]
    fn panning_with_a_selection_moves_it_one_base_without_refetching() {
        let mut app = demo_app(0);
        app.needs_fetch = false;
        let selected = app.view_start + 10;
        app.cache.loaded_region = Some(Region::new("chrDemo", 0, app.current_contig_len()));

        app.select_reference_position(selected);

        assert_eq!(app.selected_ref_pos, Some(selected));
        assert_eq!(
            app.selected_allele_tally,
            Some(PositionAlleleTally::default())
        );
        assert!(!app.needs_fetch);

        let start = app.view_start;
        app.pan(20);

        assert_eq!(app.selected_ref_pos, Some(selected + 1));
        assert_eq!(app.view_start, start);
        assert!(!app.needs_fetch);

        app.pan(-20);

        assert_eq!(app.selected_ref_pos, Some(selected));
    }

    #[test]
    fn zooming_preserves_and_centers_the_selected_base_without_refetching() {
        let mut app = demo_app(0);
        app.source.contigs[0].length = 10_000;
        app.view_start = 100;
        app.view_end = 1_100;
        app.cache.loaded_region = Some(Region::new("chrDemo", 0, 10_000));
        app.needs_fetch = false;
        app.select_reference_position(200);

        app.zoom_in();

        assert_eq!(app.selected_ref_pos, Some(200));
        assert!((app.view_start..app.view_end).contains(&200));
        assert_eq!(app.view_start, 0);
        assert!(!app.needs_fetch);

        app.zoom_out();

        assert_eq!(app.selected_ref_pos, Some(200));
        assert!((app.view_start..app.view_end).contains(&200));
        assert!(!app.needs_fetch);
    }

    #[test]
    fn selected_base_panning_activates_the_insertion_at_its_anchor() {
        let mut app = demo_app(0);
        app.view_start = 0;
        app.view_end = 20;
        app.terminal_cols = 22;
        app.expand_insertions = true;
        app.selected_insertion_ref_pos = Some(2);
        app.cache.reads = vec![
            read_with_insertion("first-insertion", 3),
            read_with_insertion("second-insertion", 10),
        ];
        app.cache.pileup_rows = vec![vec![0, 1]];
        app.selected_ref_pos = Some(8);

        app.pan(1);

        assert_eq!(app.selected_ref_pos, Some(9));
        assert_eq!(app.selected_insertion_ref_pos, Some(9));
    }

    #[test]
    fn enabling_insertion_expansion_activates_the_selected_anchor() {
        let mut app = demo_app(0);
        app.view_start = 0;
        app.view_end = 20;
        app.terminal_cols = 22;
        app.cache.reads = vec![
            read_with_insertion("first-insertion", 3),
            read_with_insertion("second-insertion", 10),
        ];
        app.cache.pileup_rows = vec![vec![0, 1]];
        app.selected_ref_pos = Some(9);

        app.toggle_insertions();

        assert!(app.expand_insertions);
        assert_eq!(app.selected_insertion_ref_pos, Some(9));
        assert_eq!(
            app.selected_insertion_gap(&app.base_view_transform())
                .map(|gap| gap.anchor_ref_pos()),
            Some(9)
        );
    }

    #[test]
    fn large_pan_is_fixed_at_one_thousand_bp() {
        let mut app = demo_app(0);
        app.source.contigs[0].length = 10_000;
        app.view_start = 100;
        app.view_end = 200;

        app.pan_large(1);

        assert_eq!(app.view_start, 1_100);
        assert_eq!(app.view_end, 1_200);
    }

    #[test]
    fn read_tracks_scroll_independently() {
        let mut app = demo_app(0);
        app.terminal_rows = 10;
        app.cache.pileup_rows = vec![Vec::new(); 7];

        app.scroll_read_track(ReadTrack::Combined, 1);

        assert_eq!(app.active_read_track, ReadTrack::Combined);
        assert_eq!(app.read_track_scroll(ReadTrack::Combined), 1);

        app.show_phasing = true;
        app.cache.phase_layout = Some(crate::cache::PhasePileupLayout {
            hp1_rows: 0..4,
            hp2_rows: 4..6,
            hp1_viewport_rows: 2,
            hp2_viewport_rows: 1,
            ..crate::cache::PhasePileupLayout::default()
        });

        app.scroll_read_track(ReadTrack::Hp2, 1);

        assert_eq!(app.active_read_track, ReadTrack::Hp2);
        assert_eq!(app.read_track_scroll(ReadTrack::Hp2), 1);
        assert_eq!(app.read_track_scroll(ReadTrack::Combined), 1);

        app.cycle_active_read_track(false);
        assert_eq!(app.active_read_track, ReadTrack::Hp1);
    }

    #[test]
    fn mapq_filter_refreshes_selected_allele_tally() {
        let mut app = demo_app(0);
        app.view_start = 0;
        app.view_end = 10;
        let mut low = crate::cache::RenderRead {
            name: "low".to_string(),
            start: 0,
            end: 3,
            strand: crate::cache::Strand::Forward,
            mapq: 20,
            cigar_ops: vec![crate::cache::CigarOp::Match(3)],
            sequence: b"AAA".to_vec(),
            methylation: Vec::new(),
            deleted_reference_sequences: Vec::new(),
            phase: Default::default(),
            is_secondary: false,
            is_supplementary: false,
            is_duplicate: false,
        };
        low.sequence[1] = b'G';
        app.cache.reads = vec![low];
        app.select_reference_position(1);
        assert_eq!(
            app.selected_allele_tally
                .as_ref()
                .and_then(|tally| tally.base_counts.get(&b'G')),
            Some(&1)
        );

        app.begin_mapq_filter();
        app.command_buffer.push_str("30");
        app.confirm_mapq_filter();

        assert_eq!(
            app.selected_allele_tally
                .as_ref()
                .map(|tally| &tally.base_counts),
            Some(&std::collections::BTreeMap::new())
        );
    }

    #[test]
    fn selected_demo_deletion_uses_read_consensus_without_a_reference() {
        let mut app = demo_app(0);
        app.terminal_cols = 6;
        app.terminal_rows = 20;
        app.jump_to_region(&Region::new("chrDemo", 61, 65))
            .expect("set deletion region");
        app.refresh().expect("load demo reads");

        app.select_reference_position(62);

        assert_eq!(
            app.selected_allele_tally
                .as_ref()
                .and_then(|tally| tally.deletion_counts.get(b"GT" as &[u8])),
            Some(&1)
        );
        assert_eq!(
            app.selected_allele_tally
                .as_ref()
                .map(|tally| tally.deletion_count),
            Some(0)
        );
    }
}
