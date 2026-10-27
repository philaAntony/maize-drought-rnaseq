library(tximport)
samples <- list.files("quants")
files <- file.path("quants", samples, "quant.sf")
names(files) <- samples
tx2gene <- read.delim("reference/tx2gene.tsv", header = FALSE,
                      col.names = c("transcript", "gene"))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
print(dim(txi$counts))
write.csv(txi$counts, "counts_matrix.csv")
