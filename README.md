# Maize drought stress RNA-seq (Project 1)

RNA-seq analysis of drought stress in maize: QC, quantification with Salmon, differential expression with DESeq2, and functional interpretation.

## Reference
- Source: Ensembl Plants release 63
- Assembly: Zm-B73-REFERENCE-NAM-5.0
- Files: cdna.all FASTA + full GTF (not chr, not abinitio)
- Check: all 72,539 FASTA transcripts are in the GTF; the 4,802 extra GTF entries are misc_non_coding RNAs, so counts cover protein-coding genes only.
