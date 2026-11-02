#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p reference

REL=63
BASE=https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-${REL}/plants
FA=Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa.gz
GTF=Zea_mays.Zm-B73-REFERENCE-NAM-5.0.${REL}.gtf.gz

[ -s reference/$FA ]  || wget -q $BASE/fasta/zea_mays/cdna/$FA -O reference/$FA
[ -s reference/$GTF ] || wget -q $BASE/gtf/zea_mays/$GTF -O reference/$GTF

# Every FASTA transcript must be in the GTF
missing=$(comm -23 \
  <(zcat reference/$FA | grep '>' | cut -d' ' -f1 | sed 's/>//' | sort -u) \
  <(zcat reference/$GTF | awk '$3=="transcript"' | grep -o 'transcript_id "[^"]*"' | cut -d'"' -f2 | sort -u) \
  | wc -l)
echo "FASTA transcripts missing from GTF: $missing (expected 0)"
[ "$missing" -eq 0 ] || { echo "FASTA and GTF do not match"; exit 1; }
