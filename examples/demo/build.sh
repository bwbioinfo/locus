#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

samtools faidx demo.fa
samtools view -bS demo.sam \
  | samtools sort -o demo.sorted.bam -
samtools index demo.sorted.bam

cargo run --quiet --manifest-path ../../Cargo.toml -- prepare-annotations demo.gff --output demo.sorted.gff.gz
