#!/usr/bin/env bash
set -euo pipefail
# Mapping rate of every run, read from Salmon's meta_info.json
mkdir -p results
tail -n +2 metadata/sample_metadata.tsv | while IFS=$'\t' read -r id gsm cond wp rep run; do
  rate=$(grep -o '"percent_mapped": [0-9.]*' quants/$run/aux_info/meta_info.json | grep -o '[0-9.]*$')
  printf "%s\t%s\t%s\n" "$id" "$run" "$rate"
done > results/mapping_rates.tsv
