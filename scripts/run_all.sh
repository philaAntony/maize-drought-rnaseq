#!/usr/bin/env bash
# Runs the whole workflow from the repository root: bash scripts/run_all.sh
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs results figures reference

step() { echo; echo "[$(date +%H:%M:%S)] $*"; }

step "00 reference download";  bash scripts/00_download_reference.sh

# Reference files must be present (downloaded above)
for pat in "reference/*.cdna.all.fa.gz" "reference/*.gtf.gz"; do
  ls $pat >/dev/null 2>&1 || { echo "Missing $pat. Download the reference files first (see README)."; exit 1; }
done

# GO annotations (MaizeGDB Pannzer predictions)
if [ ! -f reference/pannzer_B73.out ]; then
  step "Downloading Pannzer GO annotations"
  wget -q https://download.maizegdb.org/GeneFunction_and_Expression/Pannzer_GO_Terms/B73.out -O reference/pannzer_B73.out
fi

if [ "${SKIP_DOWNLOAD:-0}" = "1" ]; then step "01 download (skipped)"; else step "01 download"; bash scripts/01_download.sh; fi
step "02 raw QC";              bash scripts/02_qc_raw.sh
step "03 trimming";            bash scripts/03_trim.sh
step "04 Salmon index";        bash scripts/04_build_index.sh
step "05 quantification";      bash scripts/05_quantify.sh
step "05b mapping rates";      bash scripts/05b_mapping_rates.sh
step "06 counts matrix";       Rscript scripts/06_counts_matrix.R
step "07 DESeq2";              Rscript scripts/07_deseq2.R
step "08 volcano and heatmap"; Rscript scripts/08_volcano_heatmap.R
step "09 robustness check";    Rscript scripts/09_robustness_sev2.R
step "10a gene lists";         Rscript scripts/10a_export_gene_lists.R
step "10 GO enrichment";       Rscript scripts/10_go_enrichment.R

step "Finished"
