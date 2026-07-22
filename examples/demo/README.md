# Locus Demo Dataset

This directory contains a tiny synthetic dataset for README screenshots and manual smoke tests.

Files:

- `demo.fa`: reference sequence for `chrDemo`
- `demo.gff`: unsorted source annotation file
- `demo.sam`: source alignments with insertion, deletion, MM/ML methylation, and HP/PS phasing tags
- `build.sh`: regenerates indexed browser-ready files from the source inputs
- `capture.sh`: rebuilds the data, drives the TUI demo keys, and refreshes README screenshot artifacts

Generated browser-ready files:

- `demo.fa.fai`
- `demo.sorted.bam`
- `demo.sorted.bam.bai`
- `demo.sorted.gff.gz`
- `demo.sorted.gff.gz.tbi`

Run:

```bash
cargo run -- examples/demo/demo.sorted.bam \
  --region chrDemo:45-115 \
  --reference examples/demo/demo.fa \
  --gff examples/demo/demo.sorted.gff.gz
```

Inside the TUI, press `i`, `Tab`, `m`, `p`, and `s` to reproduce the expanded-insertion, methylation- and phasing-enabled screenshot. The fixture contains HP1, HP2, untagged, and malformed phase-tag examples.

To regenerate the committed screenshot artifacts:

```bash
examples/demo/capture.sh
```

This writes:

- `docs/captures/demo-expanded-methylation.html`
- `docs/captures/demo-expanded-methylation.ansi.txt`
- `docs/images/demo-expanded-methylation.png` when `chromium` is available
