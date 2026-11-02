
## Week 6: Quantification
- **Did:** built a Salmon index from the Ensembl cDNA, checked FASTA/GTF consistency, quantified 36 samples, built a gene-level counts matrix.
- **Broke:** the Salmon index step needed the right file pairing; mixed up two repos and pasted R code into the shell; missing R packages (jsonlite).
- **Learned:** counts vs TPM and why DESeq2 wants counts; how to verify a reference and annotation match; reproducibility means saving commands as scripts.

## Week 7: Differential expression
- **Did:** collapsed 36 sequencing runs into 9 biological samples, ran DESeq2 (mild and severe vs well-watered), made PCA, distance heatmap and MA plots.
- **Broke:** missing R packages (jsonlite); only found that the 36 SRR runs were technical replicates by reading the metadata before analysing.
- **Learned:** technical runs must be summed, not treated as replicates; PC1 (88.6%) tracked the water potential gradient; SEV_2 is the most extreme severe replicate but has normal mapping quality, so it was kept.

## Week 8 (so far): Filtered lists and robustness
- **Did:** applied |log2FC| >= 1 with padj < 0.05, made volcano plots and a heatmap, ran a leave-one-out check on SEV_2.
- **Learned:** a fold-change threshold cut the severe list from 13,729 to 4,095 genes; SEV_2 has the largest influence on the severe result, yet 78% to 82% of genes hold without it.

## Week 8: Enrichment and interpretation
- **Did:** applied a fold-change threshold, made volcano plots and a heatmap, ran a leave-one-out check on SEV_2, built a gene-to-GO table from MaizeGDB Pannzer predictions, and ran GO enrichment on up and down lists.
- **Broke:** BioMart returned HTTP 405 and a browser export never downloaded; switched to a file download from MaizeGDB.
- **Learned:** enrichment terms overlap heavily so they must be grouped into themes; redox terms turned out to be peroxidase-driven and are ambiguous without looking at the genes; "not enriched" is not the same as "not changed".

## Week 9: Reproducibility
- **Did:** restructured the README, fixed the environment file, wrote run_all.sh and helper scripts, and reran the workflow from a fresh clone.
- **Learned:** manual terminal steps must become scripts or they cannot be reproduced; a reproduction test finds problems a normal run hides; the rerun gave identical counts, DESeq2 results, gene lists and GO terms.

## Week 9: Reproducibility
- **Did:** restructured the README, fixed the environment file, wrote run_all.sh and helper scripts, and reran the workflow from a fresh clone.
- **Broke:** the first environment file lacked the bioconda channel and several R packages; the first reproduction attempt stopped silently because wget -q hid a download error.
- **Learned:** manual terminal steps must become scripts or they cannot be reproduced; a reproduction test finds problems a normal run hides; the rerun gave identical counts, DESeq2 results, gene lists and GO terms.
