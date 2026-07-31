use std::{collections::BTreeMap, ops::Range};

use crate::{reference::ReferenceSlice, region::Region};

/// Strand of a read alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

/// A CIGAR operation in our render model.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub enum CigarOp {
    Match(u64),
    Mismatch(u64),
    Insertion(u64),
    Deletion(u64),
    Skip(u64),
    SoftClip(u64),
}

impl CigarOp {
    pub fn ref_len(&self) -> u64 {
        match self {
            CigarOp::Match(n) | CigarOp::Mismatch(n) | CigarOp::Deletion(n) | CigarOp::Skip(n) => {
                *n
            }
            _ => 0,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ModificationStrand {
    Forward,
    Reverse,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ModifiedBaseCall {
    pub read_pos: usize,
    pub canonical_base: u8,
    pub strand: ModificationStrand,
    pub modification: String,
    pub probability: Option<u8>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignedModifiedBaseCall {
    pub call: ModifiedBaseCall,
    pub ref_pos: Option<u64>,
}

/// Standard haplotagging metadata carried by HP and PS alignment tags.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ReadPhase {
    pub haplotype: Option<u32>,
    pub phase_set: Option<u32>,
}

/// Lightweight read representation for rendering.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct RenderRead {
    pub name: String,
    /// 0-based start on the reference
    pub start: u64,
    /// 0-based exclusive end on the reference
    pub end: u64,
    pub strand: Strand,
    pub mapq: u8,
    pub cigar_ops: Vec<CigarOp>,
    /// ASCII-decoded read sequence (A/C/G/T/N), read-coordinate indexed.
    pub sequence: Vec<u8>,
    pub methylation: Vec<ModifiedBaseCall>,
    /// Deleted reference alleles from the MD tag, in CIGAR deletion order.
    pub deleted_reference_sequences: Vec<Vec<u8>>,
    pub phase: ReadPhase,
    pub is_secondary: bool,
    pub is_supplementary: bool,
    pub is_duplicate: bool,
}

/// Alleles observed at one reference position across a set of reads.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PositionAlleleTally {
    pub base_counts: BTreeMap<u8, usize>,
    /// Deleted reference alleles when the corresponding FASTA sequence is available.
    pub deletion_counts: BTreeMap<Vec<u8>, usize>,
    /// Deletions whose reference sequence was not available in the cached slice.
    pub deletion_count: usize,
    pub insertion_counts: BTreeMap<Vec<u8>, usize>,
}

impl PositionAlleleTally {
    fn add_base(&mut self, base: u8) {
        *self
            .base_counts
            .entry(base.to_ascii_uppercase())
            .or_default() += 1;
    }

    fn add_insertion(&mut self, sequence: Vec<u8>) {
        *self.insertion_counts.entry(sequence).or_default() += 1;
    }

    fn add_deletion(&mut self, sequence: Option<Vec<u8>>) {
        if let Some(sequence) = sequence {
            *self.deletion_counts.entry(sequence).or_default() += 1;
        } else {
            self.deletion_count += 1;
        }
    }
}

impl RenderRead {
    #[allow(dead_code)]
    pub fn len_bp(&self) -> u64 {
        self.end.saturating_sub(self.start)
    }

    pub fn aligned_methylation(&self) -> Vec<AlignedModifiedBaseCall> {
        let mut aligned = self
            .methylation
            .iter()
            .cloned()
            .map(|call| AlignedModifiedBaseCall {
                call,
                ref_pos: None,
            })
            .collect::<Vec<_>>();

        let mut read_pos: usize = 0;
        let mut ref_pos = self.start;

        for &op in &self.cigar_ops {
            match op {
                CigarOp::SoftClip(n) | CigarOp::Insertion(n) => {
                    read_pos += n as usize;
                }
                CigarOp::Match(n) | CigarOp::Mismatch(n) => {
                    let start_read_pos = read_pos;
                    let end_read_pos = read_pos + n as usize;
                    for aligned_call in &mut aligned {
                        if (start_read_pos..end_read_pos).contains(&aligned_call.call.read_pos) {
                            aligned_call.ref_pos = Some(
                                ref_pos + (aligned_call.call.read_pos - start_read_pos) as u64,
                            );
                        }
                    }
                    read_pos = end_read_pos;
                    ref_pos += n;
                }
                CigarOp::Deletion(n) | CigarOp::Skip(n) => {
                    ref_pos += n;
                }
            }
        }

        aligned
    }

    fn tally_alleles_at(
        &self,
        position: u64,
        reference: Option<&ReferenceSlice>,
        tally: &mut PositionAlleleTally,
        unresolved_deletions: &mut BTreeMap<(u64, u64), usize>,
    ) {
        let mut read_pos = 0usize;
        let mut ref_pos = self.start;
        let mut deletion_idx = 0usize;

        for &op in &self.cigar_ops {
            match op {
                CigarOp::SoftClip(n) => read_pos += n as usize,
                CigarOp::Match(n) | CigarOp::Mismatch(n) => {
                    let end = ref_pos.saturating_add(n);
                    if (ref_pos..end).contains(&position) {
                        let offset = (position - ref_pos) as usize;
                        let base = self
                            .sequence
                            .get(read_pos + offset)
                            .copied()
                            .unwrap_or(b'N');
                        tally.add_base(base);
                    }
                    read_pos += n as usize;
                    ref_pos = end;
                }
                CigarOp::Insertion(n) => {
                    if ref_pos.saturating_sub(1) == position {
                        let sequence = (0..n as usize)
                            .map(|offset| {
                                self.sequence
                                    .get(read_pos + offset)
                                    .copied()
                                    .unwrap_or(b'N')
                                    .to_ascii_uppercase()
                            })
                            .collect();
                        tally.add_insertion(sequence);
                    }
                    read_pos += n as usize;
                }
                CigarOp::Deletion(n) => {
                    let end = ref_pos.saturating_add(n);
                    if (ref_pos..end).contains(&position) {
                        let sequence = reference
                            .and_then(|reference| {
                                (ref_pos..end)
                                    .map(|deleted_pos| reference.base_at(deleted_pos))
                                    .collect::<Option<Vec<_>>>()
                            })
                            .or_else(|| {
                                self.deleted_reference_sequences
                                    .get(deletion_idx)
                                    .filter(|sequence| sequence.len() == n as usize)
                                    .cloned()
                            });
                        if let Some(sequence) = sequence {
                            tally.add_deletion(Some(sequence));
                        } else {
                            tally.add_deletion(None);
                            *unresolved_deletions.entry((ref_pos, n)).or_default() += 1;
                        }
                    }
                    deletion_idx += 1;
                    ref_pos = end;
                }
                // A CIGAR skip is an intron or other reference skip, not an indel allele.
                CigarOp::Skip(n) => ref_pos += n,
            }
        }
    }

    /// Return the query base aligned to one reference position, if this read covers it.
    fn base_at_reference_position(&self, position: u64) -> Option<u8> {
        let mut read_pos = 0usize;
        let mut ref_pos = self.start;

        for &op in &self.cigar_ops {
            match op {
                CigarOp::SoftClip(n) | CigarOp::Insertion(n) => read_pos += n as usize,
                CigarOp::Match(n) | CigarOp::Mismatch(n) => {
                    let end = ref_pos.saturating_add(n);
                    if (ref_pos..end).contains(&position) {
                        return self
                            .sequence
                            .get(read_pos + (position - ref_pos) as usize)
                            .copied();
                    }
                    read_pos += n as usize;
                    ref_pos = end;
                }
                CigarOp::Deletion(n) | CigarOp::Skip(n) => {
                    let end = ref_pos.saturating_add(n);
                    if (ref_pos..end).contains(&position) {
                        return None;
                    }
                    ref_pos = end;
                }
            }
        }

        None
    }
}

/// A single row of the pileup, containing non-overlapping reads.
/// Each entry is an index into the `reads` Vec.
pub type PileupRow = Vec<usize>;

/// Per-terminal-column coverage count.
pub type CoverageBins = Vec<u32>;

/// Row ranges for independently packed haplotype tracks.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PhasePileupLayout {
    pub hp1_rows: Range<usize>,
    pub hp2_rows: Range<usize>,
    pub unphased_rows: Option<Range<usize>>,
    /// Rows visible at once in each independently scrollable section.
    pub hp1_viewport_rows: usize,
    pub hp2_viewport_rows: usize,
    pub unphased_viewport_rows: usize,
}

/// Cached data for the currently visible region.
#[derive(Default)]
pub struct RegionCache {
    /// The padded region that was actually fetched from BAM.
    pub loaded_region: Option<Region>,
    /// All reads fetched for the padded region.
    pub reads: Vec<RenderRead>,
    /// Greedy row-packed pileup layout.
    pub pileup_rows: Vec<PileupRow>,
    /// Group boundaries when phased-read track separation is enabled.
    pub phase_layout: Option<PhasePileupLayout>,
    /// Per-column coverage histogram (length == terminal_cols).
    pub coverage: CoverageBins,
    /// Reference bases for the padded region, when a FASTA was supplied.
    pub reference: Option<ReferenceSlice>,
}

impl RegionCache {
    #[allow(dead_code)]
    pub fn clear(&mut self) {
        self.loaded_region = None;
        self.reads.clear();
        self.pileup_rows.clear();
        self.phase_layout = None;
        self.coverage.clear();
        self.reference = None;
    }

    /// Rebuild the complete pileup layout for the visible sub-region.
    ///
    /// `max_height` determines the viewport capacity; it never discards packed rows.
    pub fn layout_pileup(
        &mut self,
        visible: &Region,
        max_height: usize,
        min_mapq: u8,
        separate_phases: bool,
    ) {
        let visible_reads: Vec<usize> = self
            .reads
            .iter()
            .enumerate()
            .filter(|(_, read)| {
                read.mapq >= min_mapq && read.start < visible.end && read.end > visible.start
            })
            .map(|(i, _)| i)
            .collect();

        if separate_phases {
            self.layout_phased_pileup(&visible_reads, max_height);
        } else {
            self.pileup_rows = pack_reads(&visible_reads, &self.reads, usize::MAX);
            self.phase_layout = None;
        }
    }

    fn layout_phased_pileup(&mut self, visible_reads: &[usize], max_height: usize) {
        let mut hp1 = Vec::new();
        let mut hp2 = Vec::new();
        let mut unphased = Vec::new();

        for &read_idx in visible_reads {
            match self.reads[read_idx].phase.haplotype {
                Some(1) => hp1.push(read_idx),
                Some(2) => hp2.push(read_idx),
                _ => unphased.push(read_idx),
            }
        }

        let has_unphased = !unphased.is_empty();
        let header_rows = 2 + usize::from(has_unphased);
        let row_budget = max_height.saturating_sub(header_rows);
        let hp1_rows = pack_reads(&hp1, &self.reads, usize::MAX);
        let hp2_rows = pack_reads(&hp2, &self.reads, usize::MAX);
        let unphased_rows = pack_reads(&unphased, &self.reads, usize::MAX);
        let (hp1_viewport_rows, hp2_viewport_rows, unphased_viewport_rows) =
            reclaim_unused_phase_rows(
                row_budget,
                [hp1_rows.len(), hp2_rows.len(), unphased_rows.len()],
                has_unphased,
            );

        let hp1_end = hp1_rows.len();
        let hp2_end = hp1_end + hp2_rows.len();
        let unphased_end = hp2_end + unphased_rows.len();

        self.pileup_rows = hp1_rows;
        self.pileup_rows.extend(hp2_rows);
        self.pileup_rows.extend(unphased_rows);
        self.phase_layout = Some(PhasePileupLayout {
            hp1_rows: 0..hp1_end,
            hp2_rows: hp1_end..hp2_end,
            unphased_rows: has_unphased.then_some(hp2_end..unphased_end),
            hp1_viewport_rows,
            hp2_viewport_rows,
            unphased_viewport_rows,
        });
    }

    /// Compute per-column coverage over `visible` region, binned to `cols` columns.
    pub fn compute_coverage(&mut self, visible: &Region, cols: usize, min_mapq: u8) {
        self.coverage = bin_coverage(&self.reads, visible, cols, min_mapq);
    }

    /// Count read alleles at `position` across all loaded reads passing `min_mapq`.
    ///
    /// This intentionally does not use `pileup_rows`: row capacity is a display concern and
    /// should not change a position's allele tally.
    pub fn allele_tally_at(&self, position: u64, min_mapq: u8) -> PositionAlleleTally {
        let mut tally = PositionAlleleTally::default();
        let mut unresolved_deletions = BTreeMap::new();
        for read in self.reads.iter().filter(|read| read.mapq >= min_mapq) {
            read.tally_alleles_at(
                position,
                self.reference.as_ref(),
                &mut tally,
                &mut unresolved_deletions,
            );
        }

        for ((start, len), count) in unresolved_deletions {
            let end = start.saturating_add(len);
            let Some(sequence) = self.pileup_consensus_sequence(start, end, min_mapq) else {
                continue;
            };
            tally.deletion_count = tally.deletion_count.saturating_sub(count);
            for _ in 0..count {
                tally.add_deletion(Some(sequence.clone()));
            }
        }
        tally
    }

    /// Infer absent deleted reference bases from reads that align across the deletion span.
    fn pileup_consensus_sequence(&self, start: u64, end: u64, min_mapq: u8) -> Option<Vec<u8>> {
        (start..end)
            .map(|position| {
                let mut counts = BTreeMap::new();
                for read in self.reads.iter().filter(|read| read.mapq >= min_mapq) {
                    let Some(base) = read.base_at_reference_position(position) else {
                        continue;
                    };
                    *counts.entry(base.to_ascii_uppercase()).or_insert(0usize) += 1;
                }
                counts
                    .into_iter()
                    .max_by_key(|(_, count)| *count)
                    .map(|(base, _)| base)
            })
            .collect()
    }
}

fn phase_row_limits(total_rows: usize, has_unphased: bool) -> (usize, usize, usize) {
    if !has_unphased || total_rows < 3 {
        return (total_rows.div_ceil(2), total_rows / 2, 0);
    }

    let unphased_rows = (total_rows / 5).max(1);
    let phased_rows = total_rows - unphased_rows;
    (phased_rows.div_ceil(2), phased_rows / 2, unphased_rows)
}

fn reclaim_unused_phase_rows(
    total_rows: usize,
    required_rows: [usize; 3],
    has_unphased: bool,
) -> (usize, usize, usize) {
    let (hp1, hp2, unphased) = phase_row_limits(total_rows, has_unphased);
    let mut limits = [hp1, hp2, unphased];
    for (limit, required) in limits.iter_mut().zip(required_rows) {
        *limit = (*limit).min(required);
    }

    let mut remaining = total_rows.saturating_sub(limits.into_iter().sum::<usize>());
    while remaining > 0 {
        let mut assigned = false;
        for idx in 0..limits.len() {
            if limits[idx] >= required_rows[idx] {
                continue;
            }
            limits[idx] += 1;
            remaining -= 1;
            assigned = true;
            if remaining == 0 {
                break;
            }
        }
        if !assigned {
            break;
        }
    }

    (limits[0], limits[1], limits[2])
}

/// Greedy row-packing: sort reads by start, assign each to the first row where it fits.
fn pack_reads(indices: &[usize], reads: &[RenderRead], max_rows: usize) -> Vec<PileupRow> {
    // Sort by start coordinate
    let mut sorted = indices.to_vec();
    sorted.sort_by_key(|&i| reads[i].start);

    // row_ends[r] = exclusive end of last read placed in row r
    let mut row_ends: Vec<u64> = Vec::new();
    let mut rows: Vec<PileupRow> = Vec::new();

    for &idx in &sorted {
        let read = &reads[idx];
        // Find first row where this read doesn't overlap (with 1-col gap for readability)
        let target_row = row_ends
            .iter()
            .position(|&end| read.start > end)
            .unwrap_or(row_ends.len());

        if target_row >= max_rows {
            // The caller requested a bounded packing layout.
            continue;
        }

        if target_row == rows.len() {
            rows.push(Vec::new());
            row_ends.push(0);
        }

        rows[target_row].push(idx);
        row_ends[target_row] = read.end;
    }

    rows
}

/// Bin per-base coverage into `cols` terminal columns.
fn bin_coverage(reads: &[RenderRead], visible: &Region, cols: usize, min_mapq: u8) -> Vec<u32> {
    if cols == 0 || visible.len() == 0 {
        return vec![0; cols];
    }
    let mut bins = vec![0u32; cols];
    let region_len = visible.len() as f64;
    let bp_per_col = region_len / cols as f64;

    for read in reads {
        if read.mapq < min_mapq {
            continue;
        }

        // intersect read with visible region
        let r_start = read.start.max(visible.start);
        let r_end = read.end.min(visible.end);
        if r_start >= r_end {
            continue;
        }
        // map to columns
        let col_start = ((r_start - visible.start) as f64 / bp_per_col) as usize;
        let col_end = ((r_end - visible.start) as f64 / bp_per_col).ceil() as usize;
        let col_end = col_end.min(cols);
        for bin in bins.iter_mut().take(col_end).skip(col_start) {
            *bin = bin.saturating_add(1);
        }
    }
    bins
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_read(name: &str, start: u64, end: u64) -> RenderRead {
        let len = (end - start) as usize;
        RenderRead {
            name: name.to_string(),
            start,
            end,
            strand: Strand::Forward,
            mapq: 60,
            cigar_ops: vec![CigarOp::Match(end - start)],
            sequence: vec![b'A'; len],
            methylation: Vec::new(),
            deleted_reference_sequences: Vec::new(),
            phase: ReadPhase::default(),
            is_secondary: false,
            is_supplementary: false,
            is_duplicate: false,
        }
    }

    #[test]
    fn test_pack_reads_no_overlap() {
        let reads = vec![make_read("r1", 0, 100), make_read("r2", 200, 300)];
        let rows = pack_reads(&[0, 1], &reads, 10);
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].len(), 2);
    }

    #[test]
    fn test_pack_reads_overlap() {
        let reads = vec![make_read("r1", 0, 100), make_read("r2", 50, 150)];
        let rows = pack_reads(&[0, 1], &reads, 10);
        assert_eq!(rows.len(), 2);
    }

    #[test]
    fn test_pack_reads_max_rows() {
        let reads: Vec<RenderRead> = (0..5).map(|i| make_read("r", i * 10, i * 10 + 5)).collect();
        let indices: Vec<usize> = (0..5).collect();
        let rows = pack_reads(&indices, &reads, 2);
        // All non-overlapping, so should all fit in 1 row, but max_rows=2 is fine
        assert!(rows.len() <= 2);
    }

    #[test]
    fn test_bin_coverage_simple() {
        let visible = Region::new("chr1", 0, 100);
        let reads = vec![make_read("r1", 0, 50), make_read("r2", 50, 100)];
        let bins = bin_coverage(&reads, &visible, 10, 0);
        assert_eq!(bins.len(), 10);
        // each read covers 5 cols
        let total: u32 = bins.iter().sum();
        assert_eq!(total, 10);
    }

    #[test]
    fn test_bin_coverage_overlap() {
        let visible = Region::new("chr1", 0, 100);
        let reads = vec![make_read("r1", 0, 100), make_read("r2", 0, 100)];
        let bins = bin_coverage(&reads, &visible, 10, 0);
        assert!(bins.iter().all(|&c| c == 2));
    }

    #[test]
    fn mapq_filter_applies_to_pileup_and_coverage() {
        let visible = Region::new("chr1", 0, 100);
        let mut low = make_read("low", 0, 100);
        low.mapq = 29;
        low.phase.haplotype = Some(1);
        let mut high = make_read("high", 0, 100);
        high.mapq = 30;
        high.phase.haplotype = Some(2);
        let mut cache = RegionCache {
            reads: vec![low, high],
            ..RegionCache::default()
        };

        cache.layout_pileup(&visible, 10, 30, false);
        cache.compute_coverage(&visible, 10, 30);

        assert_eq!(cache.pileup_rows, vec![vec![1]]);
        assert!(cache.coverage.iter().all(|&count| count == 1));
        assert_eq!(cache.reads.len(), 2);

        cache.layout_pileup(&visible, 10, 30, true);
        let layout = cache.phase_layout.as_ref().expect("phase layout");
        assert!(layout.hp1_rows.is_empty());
        assert_eq!(layout.hp2_rows, 0..1);
        assert_eq!(cache.pileup_rows, vec![vec![1]]);
    }

    #[test]
    fn allele_tally_counts_bases_snvs_and_reference_indels_without_skips() {
        let mut matched = make_read("matched", 100, 103);
        matched.sequence = b"AAA".to_vec();

        let mut mismatch = make_read("snv", 100, 103);
        mismatch.cigar_ops = vec![CigarOp::Mismatch(3)];
        mismatch.sequence = b"ACA".to_vec();

        let mut inserted = make_read("inserted", 100, 103);
        inserted.cigar_ops = vec![CigarOp::Match(1), CigarOp::Insertion(2), CigarOp::Match(2)];
        inserted.sequence = b"AGGTC".to_vec();

        let mut deleted = make_read("deleted", 100, 105);
        deleted.cigar_ops = vec![CigarOp::Match(1), CigarOp::Deletion(3), CigarOp::Match(1)];
        deleted.sequence = b"AC".to_vec();

        let mut skipped = make_read("skipped", 100, 104);
        skipped.cigar_ops = vec![CigarOp::Match(1), CigarOp::Skip(2), CigarOp::Match(1)];
        skipped.sequence = b"AT".to_vec();

        let mut low_mapq = make_read("low-mapq", 100, 103);
        low_mapq.mapq = 20;
        low_mapq.sequence = b"GGG".to_vec();

        let cache = RegionCache {
            reads: vec![matched, mismatch, inserted, deleted, skipped, low_mapq],
            reference: Some(ReferenceSlice {
                start: 100,
                bases: b"AACTG".to_vec(),
            }),
            ..RegionCache::default()
        };

        let tally = cache.allele_tally_at(101, 30);

        assert_eq!(tally.base_counts.get(&b'A'), Some(&1));
        assert_eq!(tally.base_counts.get(&b'C'), Some(&1));
        assert_eq!(tally.base_counts.get(&b'T'), Some(&1));
        assert!(!tally.base_counts.contains_key(&b'G'));
        assert_eq!(tally.deletion_counts.get(b"ACT" as &[u8]), Some(&1));
        assert_eq!(tally.deletion_count, 0);
        assert!(tally.insertion_counts.is_empty());

        let insertion_tally = cache.allele_tally_at(100, 30);
        assert_eq!(
            insertion_tally.insertion_counts.get(b"GG" as &[u8]),
            Some(&1)
        );
    }

    #[test]
    fn allele_tally_uses_md_deleted_sequence_without_a_reference_slice() {
        let mut deleted = make_read("deleted", 100, 106);
        deleted.cigar_ops = vec![CigarOp::Match(1), CigarOp::Deletion(4), CigarOp::Match(1)];
        deleted.sequence = b"AC".to_vec();
        deleted.deleted_reference_sequences = vec![b"ATCG".to_vec()];

        let cache = RegionCache {
            reads: vec![deleted],
            ..RegionCache::default()
        };

        let tally = cache.allele_tally_at(101, 0);

        assert_eq!(tally.deletion_counts.get(b"ATCG" as &[u8]), Some(&1));
        assert_eq!(tally.deletion_count, 0);
    }

    #[test]
    fn allele_tally_uses_pileup_consensus_for_deleted_sequence_without_reference_or_md() {
        let mut deleted = make_read("deleted", 100, 106);
        deleted.cigar_ops = vec![CigarOp::Match(1), CigarOp::Deletion(4), CigarOp::Match(1)];
        deleted.sequence = b"AC".to_vec();

        let mut support_a = make_read("support-a", 100, 106);
        support_a.sequence = b"AATCGC".to_vec();
        let mut support_b = make_read("support-b", 100, 106);
        support_b.sequence = b"AATCGC".to_vec();

        let cache = RegionCache {
            reads: vec![deleted, support_a, support_b],
            ..RegionCache::default()
        };

        let tally = cache.allele_tally_at(101, 0);

        assert_eq!(tally.deletion_counts.get(b"ATCG" as &[u8]), Some(&1));
        assert_eq!(tally.deletion_count, 0);
    }

    #[test]
    fn allele_tally_keeps_generic_deletion_when_no_sequence_source_is_available() {
        let mut deleted = make_read("deleted", 100, 106);
        deleted.cigar_ops = vec![CigarOp::Match(1), CigarOp::Deletion(4), CigarOp::Match(1)];
        deleted.sequence = b"AC".to_vec();

        let cache = RegionCache {
            reads: vec![deleted],
            ..RegionCache::default()
        };

        let tally = cache.allele_tally_at(101, 0);

        assert!(tally.deletion_counts.is_empty());
        assert_eq!(tally.deletion_count, 1);
    }

    #[test]
    fn phased_pileup_packs_haplotypes_into_independent_row_ranges() {
        let visible = Region::new("chr1", 0, 100);
        let mut hp1_a = make_read("hp1-a", 0, 60);
        hp1_a.phase.haplotype = Some(1);
        let mut hp1_b = make_read("hp1-b", 20, 80);
        hp1_b.phase.haplotype = Some(1);
        let mut hp2_a = make_read("hp2-a", 0, 60);
        hp2_a.phase.haplotype = Some(2);
        let mut hp2_b = make_read("hp2-b", 20, 80);
        hp2_b.phase.haplotype = Some(2);
        let unphased_a = make_read("unphased-a", 0, 10);
        let unphased_b = make_read("unphased-b", 20, 30);
        let mut cache = RegionCache {
            reads: vec![hp1_a, hp1_b, hp2_a, hp2_b, unphased_a, unphased_b],
            ..RegionCache::default()
        };

        cache.layout_pileup(&visible, 9, 0, true);

        let layout = cache.phase_layout.as_ref().expect("phase layout");
        assert_eq!(layout.hp1_rows, 0..2);
        assert_eq!(layout.hp2_rows, 2..4);
        assert_eq!(layout.unphased_rows, Some(4..5));
        assert_eq!(cache.pileup_rows[0], vec![0]);
        assert_eq!(cache.pileup_rows[1], vec![1]);
        assert_eq!(cache.pileup_rows[2], vec![2]);
        assert_eq!(cache.pileup_rows[3], vec![3]);
        assert_eq!(cache.pileup_rows[4], vec![4, 5]);
    }

    #[test]
    fn phased_pileup_keeps_rows_beyond_each_section_viewport() {
        let visible = Region::new("chr1", 0, 100);
        let mut hp1_a = make_read("hp1-a", 0, 100);
        hp1_a.phase.haplotype = Some(1);
        let mut hp1_b = make_read("hp1-b", 0, 100);
        hp1_b.phase.haplotype = Some(1);
        let mut hp2 = make_read("hp2", 0, 100);
        hp2.phase.haplotype = Some(2);
        let unphased = make_read("unphased", 0, 100);
        let mut cache = RegionCache {
            reads: vec![hp1_a, hp1_b, hp2, unphased],
            ..RegionCache::default()
        };

        cache.layout_pileup(&visible, 5, 0, true);

        let layout = cache.phase_layout.as_ref().expect("phase layout");
        assert_eq!(layout.hp1_rows, 0..2);
        assert_eq!(layout.hp2_rows, 2..3);
        assert_eq!(layout.unphased_rows, Some(3..4));
        assert_eq!(layout.hp1_viewport_rows, 1);
        assert_eq!(layout.hp2_viewport_rows, 1);
        assert_eq!(layout.unphased_viewport_rows, 0);
    }

    #[test]
    fn phased_pileup_reclaims_spare_rows_for_viewports() {
        let visible = Region::new("chr1", 0, 100);
        let mut hp1_a = make_read("hp1-a", 0, 100);
        hp1_a.phase.haplotype = Some(1);
        let mut hp1_b = make_read("hp1-b", 0, 100);
        hp1_b.phase.haplotype = Some(1);
        let mut hp1_c = make_read("hp1-c", 0, 100);
        hp1_c.phase.haplotype = Some(1);
        let mut hp2 = make_read("hp2", 0, 100);
        hp2.phase.haplotype = Some(2);
        let unphased = make_read("unphased", 0, 100);
        let mut cache = RegionCache {
            reads: vec![hp1_a, hp1_b, hp1_c, hp2, unphased],
            ..RegionCache::default()
        };

        cache.layout_pileup(&visible, 8, 0, true);

        let layout = cache.phase_layout.as_ref().expect("phase layout");
        assert_eq!(layout.hp1_rows, 0..3);
        assert_eq!(layout.hp2_rows, 3..4);
        assert_eq!(layout.unphased_rows, Some(4..5));
        assert_eq!(layout.hp1_viewport_rows, 3);
        assert_eq!(layout.hp2_viewport_rows, 1);
        assert_eq!(layout.unphased_viewport_rows, 1);
    }

    #[test]
    fn standard_pileup_keeps_rows_beyond_the_viewport() {
        let visible = Region::new("chr1", 0, 100);
        let cache_reads = vec![make_read("first", 0, 100), make_read("second", 0, 100)];
        let mut cache = RegionCache {
            reads: cache_reads,
            ..RegionCache::default()
        };

        cache.layout_pileup(&visible, 1, 0, false);

        assert_eq!(cache.pileup_rows, vec![vec![0], vec![1]]);
    }

    #[test]
    fn phase_row_limits_prioritize_both_haplotypes_in_tight_views() {
        assert_eq!(phase_row_limits(0, true), (0, 0, 0));
        assert_eq!(phase_row_limits(1, true), (1, 0, 0));
        assert_eq!(phase_row_limits(2, true), (1, 1, 0));
        assert_eq!(phase_row_limits(6, true), (3, 2, 1));
        assert_eq!(phase_row_limits(6, false), (3, 3, 0));
    }

    #[test]
    fn reclaim_unused_phase_rows_uses_every_available_row() {
        assert_eq!(reclaim_unused_phase_rows(5, [3, 1, 1], true), (3, 1, 1));
        assert_eq!(reclaim_unused_phase_rows(6, [4, 1, 1], true), (4, 1, 1));
    }

    fn methylated_call(read_pos: usize) -> ModifiedBaseCall {
        ModifiedBaseCall {
            read_pos,
            canonical_base: b'C',
            strand: ModificationStrand::Forward,
            modification: "m".to_string(),
            probability: Some(200),
        }
    }

    #[test]
    fn aligned_methylation_maps_match_positions() {
        let mut read = make_read("r", 100, 105);
        read.sequence = b"ACGTC".to_vec();
        read.methylation = vec![methylated_call(1), methylated_call(4)];

        let aligned = read.aligned_methylation();

        assert_eq!(aligned[0].ref_pos, Some(101));
        assert_eq!(aligned[1].ref_pos, Some(104));
    }

    #[test]
    fn aligned_methylation_respects_indels_skips_and_soft_clips() {
        let mut read = make_read("r", 100, 108);
        read.cigar_ops = vec![
            CigarOp::SoftClip(1),
            CigarOp::Match(2),
            CigarOp::Insertion(1),
            CigarOp::Match(2),
            CigarOp::Deletion(2),
            CigarOp::Skip(1),
            CigarOp::Match(1),
        ];
        read.sequence = b"SACGTA".to_vec();
        read.methylation = vec![
            methylated_call(0),
            methylated_call(2),
            methylated_call(3),
            methylated_call(4),
            methylated_call(5),
        ];

        let aligned = read.aligned_methylation();
        let ref_positions = aligned.iter().map(|call| call.ref_pos).collect::<Vec<_>>();

        assert_eq!(
            ref_positions,
            vec![None, Some(101), None, Some(102), Some(103)]
        );
    }
}
