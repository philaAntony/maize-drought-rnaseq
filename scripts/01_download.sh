#!/usr/bin/env bash
# Download runs from SRA and convert to gzipped FASTQ (single-end).
# Retries prefetch on network failures. Safe to re-run: finished runs are skipped.
# Usage: bash scripts/01_download.sh [SAMPLE_ID]   (e.g. WW_1; omit for all samples)
set -euo pipefail
SAMPLE="${1:-}"; THREADS=8; META=metadata/sample_metadata.tsv
mkdir -p data/sra data/raw

if [ -n "$SAMPLE" ]; then
  RUNS=$(awk -F'\t' -v s="$SAMPLE" 'NR>1 && $1==s {print $6}' "$META")
else
  RUNS=$(awk -F'\t' 'NR>1 {print $6}' "$META")
fi

for run in $RUNS; do
  if [ -s "data/raw/${run}.fastq.gz" ]; then echo "skip $run (already done)"; continue; fi
  echo "== $run =="
  ok=0
  for attempt in 1 2 3 4 5; do
    if prefetch -O data/sra "$run"; then ok=1; break; fi
    echo "prefetch failed for $run (attempt $attempt of 5); retrying in 60 s"
    sleep 60
  done
  if [ "$ok" -ne 1 ]; then echo "GIVING UP on $run after 5 attempts"; exit 1; fi
  fasterq-dump --threads "$THREADS" --temp data/sra -O data/raw "data/sra/$run/$run.sra"
  pigz -p "$THREADS" "data/raw/${run}.fastq"
  rm -rf "data/sra/$run"
done
echo "All done."
