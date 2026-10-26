#!/usr/bin/env bash
# FastQC on raw reads. Usage: bash scripts/02_qc_raw.sh
set -euo pipefail
mkdir -p qc/raw
fastqc -t 8 --outdir qc/raw data/raw/*.fastq.gz
