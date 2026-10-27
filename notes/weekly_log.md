
## Week 6: Quantification
- **Did:** built a Salmon index from the Ensembl cDNA, checked FASTA/GTF consistency, quantified 36 samples, built a gene-level counts matrix.
- **Broke:** the Salmon index step needed the right file pairing; mixed up two repos and pasted R code into the shell; missing R packages (jsonlite).
- **Learned:** counts vs TPM and why DESeq2 wants counts; how to verify a reference and annotation match; reproducibility means saving commands as scripts.

## Week 7: Differential expression
- **Did:** collapsed 36 sequencing runs into 9 biological samples, ran DESeq2 (mild and severe vs well-watered), made PCA, distance heatmap and MA plots.
- **Broke:** missing R packages (jsonlite); only found that the 36 SRR runs were technical replicates by reading the metadata before analysing.
- **Learned:** technical runs must be summed, not treated as replicates; PC1 (88.6%) tracked the water potential gradient; SEV_2 is the most extreme severe replicate but has normal mapping quality, so it was kept.
