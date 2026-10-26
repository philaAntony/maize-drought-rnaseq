# Project 1: Maize primary root growth zone under water deficit

## Question
Which genes and pathways change in the B73 primary root growth zone under mild (-0.3 MPa) and severe (-1.6 MPa) water deficit, compared with well-watered controls (-0.02 MPa), and do the two stress levels trigger different responses?

## Data
- GEO series: GSE128896 (BioProject PRJNA529231, SRA study SRP189504)
- Organism and line: Zea mays, B73 inbred
- Tissue: 1-cm primary root growth zone of seedlings, sampled 26 h after treatment began
- Design: 3 conditions x 3 biological replicates = 9 samples
- Sequencing: NextSeq 500, single-end 75 bp; each sample was sequenced in 4 runs (36 runs total), which are merged per sample at quantification
- Reference: B73 transcriptome from Ensembl Plants (record the annotation version used)
- Sample metadata: metadata/sample_metadata.tsv. Replicate numbers (1-3) were assigned in GSM order within each condition, as the SRA table does not give them.

## Comparisons (DESeq2)
- Mild vs well-watered
- Severe vs well-watered
- Optional: severe vs mild

## Expected outputs
- QC: MultiQC report and trimming summary
- Quantification: Salmon counts matrix
- Replicate checks: PCA and sample distance heatmap
- Differential expression: MA and volcano plots, top-gene heatmaps
- Interpretation: GO/KEGG enrichment and a short biological interpretation

## Risks and notes
- Sample labels in SRA need careful matching to treatments
- Pooled samples may cluster differently from independent replicates; check in the PCA
- Mild stress may give few differentially expressed genes, which is itself a result
- The GEO record title gives -0.03 MPa for well-watered but the abstract gives -0.02 MPa; -0.02 is used here
- The linked paper still needs confirming from the GEO record's citation field
