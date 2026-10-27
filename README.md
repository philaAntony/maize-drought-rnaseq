# Maize drought stress RNA-seq (Project 1)

RNA-seq analysis of drought stress in maize: QC, quantification with Salmon, differential expression with DESeq2, and functional interpretation.

## Reference
- Source: Ensembl Plants release 63
- Assembly: Zm-B73-REFERENCE-NAM-5.0
- Files: cdna.all FASTA + full GTF (not chr, not abinitio)
- Check: all 72,539 FASTA transcripts are in the GTF; the 4,802 extra GTF entries are misc_non_coding RNAs, so counts cover protein-coding genes only.

## Quantification (Week 6)
- **Tool:** Salmon, quasi-mapping against the Ensembl Plants release 63 cDNA (k = 31; reads are 75 bp, single-end).
- **Library type:** unstranded (`U`) in all 36 samples, detected automatically with `-l A`.
- **Mapping rate:** about 90% [add the range across all 36 samples once you've checked].
- **Output:** gene-level counts matrix (`counts_matrix.csv`), built with tximport using `reference/tx2gene.tsv`. [add dimensions, e.g. N genes x 36 samples]
- **Scripts:** `04_build_index.sh`, `05_quantify.sh`, `06_counts_matrix.R`.

### Counts vs TPM
DESeq2 takes raw (estimated) counts rather than TPM. It models counts with a negative binomial distribution, so the amount of noise depends on how many reads support each measurement, and it applies its own size-factor normalisation for sequencing depth. TPM has already been rescaled to sum to a million per sample, which throws away the read-depth information that DESeq2 relies on. TPM is useful for plotting and for comparing genes within a sample.

### Limitations
- Counts cover protein-coding genes only, because the cDNA FASTA excludes non-coding RNAs (4,802 `misc_non_coding` transcripts are in the GTF but not the FASTA).
- The reference is the B73 NAM-5.0 assembly, which may differ from the genotype used in the study.
