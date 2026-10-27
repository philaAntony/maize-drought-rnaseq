#!/usr/bin/env bash
set -euo pipefail

# Reference: Ensembl Plants release 63, Zea mays Zm-B73-REFERENCE-NAM-5.0
REF=reference/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa.gz
GTF=reference/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.63.gtf.gz

# Build Salmon index (use -k 21 if reads are shorter than 75 bp)
salmon index -t "$REF" -i salmon_index -k 31

# Build transcript-to-gene table for tximport (Week 7)
zcat "$GTF" | awk '$3=="transcript"' \
  | sed -E 's/.*gene_id "([^"]+)".*transcript_id "([^"]+)".*/\2\t\1/' \
  | sort -u > reference/tx2gene.tsv
