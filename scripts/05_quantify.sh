#!/usr/bin/env bash
set -euo pipefail

# Single-end reads, trimmed with fastp (scripts/03_trim.sh)
TRIMDIR=data/trimmed
mkdir -p quants logs

for fq in "$TRIMDIR"/*.trim.fastq.gz; do
  s=$(basename "$fq" .trim.fastq.gz)
  if [ -f quants/"$s"/quant.sf ]; then
    echo "Skipping $s (already done)"
    continue
  fi
  echo "Quantifying $s"
  salmon quant -i salmon_index -l A \
    -r "$fq" -p 4 --validateMappings \
    -o quants/"$s" 2> logs/salmon_"$s".log
done
