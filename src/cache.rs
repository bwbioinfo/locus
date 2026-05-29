use crate::region::Region;

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
    pub is_secondary: bool,
    pub is_supplementary: bool,
    pub is_duplicate: bool,
}

impl RenderRead {
    #[allow(dead_code)]
    pub fn len_bp(&self) -> u64 {
        self.end.saturating_sub(self.start)
    }
}

/// A single row of the pileup, containing non-overlapping reads.
/// Each entry is an index into the `reads` Vec.
pub type PileupRow = Vec<usize>;

/// Per-terminal-column coverage count.
pub type CoverageBins = Vec<u32>;

/// Cached data for the currently visible region.
#[derive(Default)]
pub struct RegionCache {
    /// The padded region that was actually fetched from BAM.
    pub loaded_region: Option<Region>,
    /// All reads fetched for the padded region.
    pub reads: Vec<RenderRead>,
    /// Greedy row-packed pileup layout.
    pub pileup_rows: Vec<PileupRow>,
    /// Per-column coverage histogram (length == terminal_cols).
    pub coverage: CoverageBins,
    /// How many reads were hidden because pileup depth was exceeded.
    pub hidden_reads: usize,
}

impl RegionCache {
    #[allow(dead_code)]
    pub fn clear(&mut self) {
        self.loaded_region = None;
        self.reads.clear();
        self.pileup_rows.clear();
        self.coverage.clear();
        self.hidden_reads = 0;
    }

    /// Rebuild pileup layout for the visible sub-region, limited to `max_rows`.
    pub fn layout_pileup(&mut self, visible: &Region, max_rows: usize) {
        // Filter reads to those overlapping the visible region
        let visible_reads: Vec<usize> = self
            .reads
            .iter()
            .enumerate()
            .filter(|(_, r)| r.start < visible.end && r.end > visible.start)
            .map(|(i, _)| i)
            .collect();

        self.pileup_rows = pack_reads(&visible_reads, &self.reads, max_rows);
        self.hidden_reads =
            visible_reads.len() - self.pileup_rows.iter().map(|row| row.len()).sum::<usize>();
    }

    /// Compute per-column coverage over `visible` region, binned to `cols` columns.
    pub fn compute_coverage(&mut self, visible: &Region, cols: usize) {
        self.coverage = bin_coverage(&self.reads, visible, cols);
    }
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
            .position(|&end| read.start >= end + 1)
            .unwrap_or(row_ends.len());

        if target_row >= max_rows {
            // Skip — hidden reads counted by caller
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
fn bin_coverage(reads: &[RenderRead], visible: &Region, cols: usize) -> Vec<u32> {
    if cols == 0 || visible.len() == 0 {
        return vec![0; cols];
    }
    let mut bins = vec![0u32; cols];
    let region_len = visible.len() as f64;
    let bp_per_col = region_len / cols as f64;

    for read in reads {
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
        for c in col_start..col_end {
            bins[c] = bins[c].saturating_add(1);
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
        let bins = bin_coverage(&reads, &visible, 10);
        assert_eq!(bins.len(), 10);
        // each read covers 5 cols
        let total: u32 = bins.iter().sum();
        assert_eq!(total, 10);
    }

    #[test]
    fn test_bin_coverage_overlap() {
        let visible = Region::new("chr1", 0, 100);
        let reads = vec![make_read("r1", 0, 100), make_read("r2", 0, 100)];
        let bins = bin_coverage(&reads, &visible, 10);
        assert!(bins.iter().all(|&c| c == 2));
    }
}
