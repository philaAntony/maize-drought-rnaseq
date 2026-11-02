# Maize drought stress RNA-seq (Project 1)

RNA-seq analysis of drought stress in maize: QC, quantification with Salmon, differential expression with DESeq2, and functional interpretation.

**Status:** QC, trimming, quantification and differential expression (DESeq2) are complete. Functional interpretation (volcano plots, GO/KEGG enrichment) is in progress.

## Reference
- Source: Ensembl Plants release 63
- Assembly: Zm-B73-REFERENCE-NAM-5.0
- Files: cdna.all FASTA + full GTF (not chr, not abinitio)
- Check: all 72,539 FASTA transcripts are in the GTF; the 4,802 extra GTF entries are misc_non_coding RNAs, so counts cover protein-coding genes only.

## Quantification (Week 6)
- **Tool:** Salmon, selective alignment (`--validateMappings`) against the Ensembl Plants release 63 cDNA (k = 31; reads are 75 bp, single-end).
- **Library type:** unstranded (`U`) in all 36 runs, detected automatically with `-l A`.
- **Mapping rate:** 90.2% to 91.5% across all 36 runs.
- **Output:** gene-level counts matrix (`counts_matrix.csv`), built with tximport using `reference/tx2gene.tsv`. 39091 genes x 36 runs (summed to 9 samples in Week 7).
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

### Limitations of the differential expression analysis
- With 3 replicates and half the transcriptome changing in the severe comparison, a large number of genes reach significance, so a fold-change threshold is applied in the interpretation step.
- DESeq2 assumes most genes do not change, which is strained when a large fraction of genes respond.

## Data
- **Study:** [Water-deficit responsive microRNAs in the primary root growth zone of Maize, BMC Plant Biology, from the GEO page for GSM3688038]
- **Accession:** GEO samples GSM3688038 to GSM3688046; SRA runs SRR8788496 to SRR8788531.
- **Design:** maize (B73 reference) under three water regimes: well-watered (-0.02 MPa), mild water deficit (-0.3 MPa) and severe water deficit (-1.6 MPa), with 3 biological replicates each.
- **Sequencing:** 75 bp single-end reads; each biological sample was sequenced in 4 runs (36 runs in total).

### Filtered gene lists (Week 8)
Genes were called differentially expressed at adjusted p < 0.05 and |log2 fold change| >= 1 (`scripts/08_volcano_heatmap.R`):
- Mild vs well-watered: 506 up, 289 down.
- Severe vs well-watered: 3,402 up, 693 down.
Fold changes are not shrunken, so values for low-count genes may be inflated. Volcano plots and a heatmap of the top 50 genes are in `figures/`.

### Robustness check (SEV_2)
SEV_2 was the most extreme severe replicate on the PCA. To check that it does not drive the severe vs well-watered result, the model was refit with each severe replicate removed in turn (`scripts/09_robustness_sev2.R`, results in `results/robustness_leave_one_out.csv`). Removing SEV_2 retained 78.1% of up-regulated and 82.0% of down-regulated genes, compared with 93.0% and 91.2% when SEV_1 was removed and 91.4% and 90.9% when SEV_3 was removed. Fold changes stayed highly correlated in all three refits (Spearman 0.974 to 0.988). SEV_2 therefore contributes more to the severe signal than the other replicates, but it does not drive the result, and it was retained. Each refit has only two severe replicates, so some loss of power is expected.

### Functional enrichment (Week 8)
GO enrichment (`scripts/10_go_enrichment.R`) was run separately for up- and down-regulated genes in each comparison against the 26,120 tested genes as background. Annotations are MaizeGDB Pannzer ARGOT predictions for B73 v5 (PPV >= 0.5), with parent terms added; 19,765 of the tested genes (75.7%) had at least one term. Significant terms: adjusted p < 0.05.

- **Severe water deficit, up-regulated:** responses to water, osmotic stress, abscisic acid and salt; phenylpropanoid, lignin, terpenoid and isoprenoid metabolism; amino acid and carbohydrate transport.
- **Severe water deficit, down-regulated:** hydrogen peroxide and ROS metabolism, detoxification, nitrate transport, xyloglucan metabolism and cell wall organisation.
- **Mild water deficit:** a partly overlapping but smaller response, with glutathione and sulfur metabolism and xenobiotic export up, and peroxide, detoxification and nitrate terms down.

Many GO terms share the same genes, so related terms were interpreted as themes rather than counted as independent findings.

### Limitations of the enrichment analysis
- GO annotations are computational predictions, so enriched terms are hypotheses rather than evidence of mechanism.
- 24.3% of tested genes had no annotation at the chosen confidence cutoff.
- Adjusted p-values are corrected within each gene list and ontology, not across all analyses.

### Interpretation (Week 8)
- **Water and ABA signalling:** genes up-regulated under severe stress are enriched for responses to water, osmotic stress and abscisic acid, which fits the known role of ABA in drought. Salt and osmotic terms also appear in the mild comparison, so the osmotic response starts early.
- **Secondary metabolism and the cell wall:** phenylpropanoid, lignin and terpenoid metabolism are up, while xyloglucan metabolism and cell wall organisation are down in severe stress. This is consistent with a shift from growth-related wall remodelling towards stiffer, lignified walls, though the enrichment does not show this directly.
- **Peroxidase-driven redox terms:** the down-regulated "hydrogen peroxide", "oxidative stress" and "detoxification" terms are driven by peroxidase genes. Many maize peroxidases are class III, which act in the cell wall as well as in antioxidant defence, so these terms should not be read as reduced antioxidant capacity. In the severe comparison, 22 genes with a peroxidase annotation were up-regulated and 20 were down-regulated (counted by gene description, so this mixes peroxidase types), so members of this family respond in both directions.
- **Nitrate transport** terms are down in both comparisons, and amino acid and carbohydrate transport terms are up in severe stress.
- **Photosynthesis-related terms were not enriched**, as expected for root tissue.

These are hypotheses from predicted annotations. The samples are root growth zones at 26 hours, so the results describe the early response of the root tip and should not be extended to leaves or to longer-term acclimation.
