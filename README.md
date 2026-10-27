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
- **Mapping rate:** 90.2% to 91.5% across all 36 runs.
- **Output:** gene-level counts matrix (`counts_matrix.csv`), built with tximport using `reference/tx2gene.tsv`. NGENES genes x 36 samples.
- **Scripts:** `04_build_index.sh`, `05_quantify.sh`, `06_counts_matrix.R`.

### Counts vs TPM
DESeq2 takes raw (estimated) counts rather than TPM. It models counts with a negative binomial distribution, so the amount of noise depends on how many reads support each measurement, and it applies its own size-factor normalisation for sequencing depth. TPM has already been rescaled to sum to a million per sample, which throws away the read-depth information that DESeq2 relies on. TPM is useful for plotting and for comparing genes within a sample.

### Limitations
- Counts cover protein-coding genes only, because the cDNA FASTA excludes non-coding RNAs (4,802 `misc_non_coding` transcripts are in the GTF but not the FASTA).
- The reference is the B73 NAM-5.0 assembly, which may differ from the genotype used in the study.

## Differential expression (Week 7)
- **Design:** DESeq2 with `~ condition` on 9 biological samples (3 conditions x 3 replicates). The 4 sequencing runs per sample were summed before analysis, so technical runs are not treated as replicates.
- **Conditions:** well-watered (-0.02 MPa, reference), mild water deficit (-0.3 MPa), severe water deficit (-1.6 MPa).
- **Genes tested:** 26,137 after filtering (at least 10 reads in at least 2 samples).
- **Results (adjusted p < 0.05):**
  - Mild vs well-watered: 1,442 up, 1,213 down.
  - Severe vs well-watered: 7,035 up, 6,694 down.
- **Script:** `scripts/07_deseq2.R`. Results tables are in `results/`, figures in `figures/`.

### Quality control
PCA of variance-stabilised counts (top 500 genes) separates the three conditions along PC1 (88.6% of variance), ordered by water potential: well-watered, mild, severe. PC2 explains 7.7%. Replicates cluster within conditions. One severe replicate (SEV_2) lies further along PC1 than the other two, with normal mapping quality (90.3%), so it was retained.

### Limitations
- With 3 replicates and half the transcriptome changing in the severe comparison, a large number of genes reach significance, so a fold-change threshold is applied in the interpretation step.
- DESeq2 assumes most genes do not change, which is strained when a large fraction of genes respond.
