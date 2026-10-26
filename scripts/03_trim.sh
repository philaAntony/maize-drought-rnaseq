#!/usr/bin/env bash
# Trim single-end reads with fastp, then run FastQC on the trimmed reads.
# Usage: bash scripts/03_trim.sh
set -euo pipefail
mkdir -p data/trimmed qc/fastp qc/trimmed
for f in data/raw/*.fastq.gz; do
  run=$(basename "$f" .fastq.gz)
  fastp -i "$f" -o "data/trimmed/${run}.trim.fastq.gz" \
        --qualified_quality_phred 20 --length_required 36 --thread 4 \
        -j "qc/fastp/${run}.fastp.json" -h "qc/fastp/${run}.fastp.html"
done
fastqc -t 8 --outdir qc/trimmed data/trimmed/*.fastq.gz
