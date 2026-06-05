# locus

A fast terminal genome browser for BAM files.

## Install

```bash
cargo build --release
# binary at target/release/locus
```

Requires Rust 1.85+ (edition 2024).

## Usage

```bash
# Open a BAM (auto-detects sample.bam.bai)
locus sample.bam

# Jump to a region on startup
locus sample.bam --region chr1:100000-101000
locus sample.bam --region "chr1:1,000,000-1,001,000"

# With a reference (reserved for future mismatch highlighting)
locus sample.bam --reference hg38.fa

# With annotations for feature rendering and gene search
locus sample.bam --gff hg38.ncbiRefSeq.gtf.gz
```

The BAM must be coordinate-sorted and indexed (`.bai` file beside it).
Annotation files can be GFF3 or GTF, plain text or gzip/BGZF-compressed.

## Keybindings

| Key | Action |
|-----|--------|
| `q` | Quit |
| `h` / `←` | Pan left (small step) |
| `l` / `→` | Pan right (small step) |
| `H` | Pan left (large step) |
| `L` | Pan right (large step) |
| `+` / `=` | Zoom in |
| `-` | Zoom out |
| `g` | Go to region (e.g. `chr1:100000-200000`) |
| `c` | Contig selector |
| `r` | Refresh current region |
| `s` | Save ANSI text and HTML screenshots to `screenshots/` |
| `?` | Toggle help overlay |
| `Esc` | Cancel input |
| `Enter` | Confirm input |

## UI Layout

```
┌─────────────────────────────────────────────────────────┐
│ locus  sample.bam  chr1  100000-101000   12.5 bp/col  42 reads │  ← top bar
├─────────────────────────────────────────────────────────┤
│ 100,000        100,500        101,000                   │  ← coordinate ruler
├─────────────────────────────────────────────────────────┤
│ ▄▄▅▅▆▆█████▆▆▅▅▄▄▃▃▂▂▁▁                               │  ← coverage histogram
├─────────────────────────────────────────────────────────┤
│ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>--<<<<<<<<<              │  ← read pileup
│ >>>>>>>>>>>>>>>>>>>>>>>>>>                              │
│ >>>>>>>>>>>>>>>>>>>X>>>>>>>>>>>>>>>>>>>>>>>>>           │
└─────────────────────────────────────────────────────────┘
│ q:quit  h/l:pan  H/L:big pan  +/-:zoom  g:goto  ?:help │  ← bottom bar
```

## Read Rendering

Reads are colored by mapping quality:
- **Green**: MAPQ ≥ 60
- **Light green**: MAPQ ≥ 30
- **Yellow**: MAPQ ≥ 10
- **Gray**: MAPQ < 10

CIGAR operations:
- `>` / `<` — alignment match (forward / reverse strand)
- `X` — sequence mismatch
- `I` — insertion into reference
- `-` — deletion from reference
- `~` — skip / intron (N)
- `S` — soft clip

Feature track:
- `─>─` / `─<─` — transcript or gene backbone, including intronic span
- `█` — exon
- `▓` — CDS
- `▒` — UTR

## Architecture

```
src/
├── main.rs          Entry point, terminal setup, main loop
├── cli.rs           Clap argument parsing
├── app.rs           App state, navigation logic
├── bam.rs           BamSource: open BAM+index, fetch reads
├── cache.rs         RenderRead, RegionCache, pileup layout, coverage binning
├── region.rs        Region type, region string parser
├── events.rs        Keyboard event dispatch
├── gff.rs           GFF3/GTF feature loading, parsing, and search
├── screenshot.rs    ANSI text and HTML screenshot export
├── ui.rs            ratatui frame drawing (top bar, overlays)
├── error.rs         LocusError enum
└── render/
    ├── mod.rs       ViewTransform (bp ↔ column mapping)
    ├── ruler.rs     Coordinate ruler widget
    ├── coverage.rs  Coverage histogram widget
    └── reads.rs     Read pileup widget
```

### Extending with new tracks

Implement the `Track` trait (future):

```rust
trait Track {
    fn name(&self) -> &str;
    fn load_region(&mut self, region: &Region) -> Result<()>;
    fn render(&self, frame: &mut Frame, area: Rect, view: &ViewTransform);
}
```

## Non-goals (first pass)

CRAM, VCF, BED tracks, split-read/SV visualization, methylation,
mouse interaction, remote files, multi-sample browsing.
