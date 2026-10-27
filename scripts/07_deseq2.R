library(tximport)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# ---- Metadata (one row per run) ----
runs <- read.delim("metadata/sample_metadata.tsv", stringsAsFactors = FALSE)
files <- file.path("quants", runs$run, "quant.sf")
names(files) <- runs$run
stopifnot(all(file.exists(files)))

# ---- Import at gene level (per run) ----
tx2gene <- read.delim("reference/tx2gene.tsv", header = FALSE,
                      col.names = c("transcript", "gene"))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# ---- Collapse technical runs: sum counts, average lengths/abundance ----
grp <- factor(runs$sample_id, levels = unique(runs$sample_id))
n_runs <- table(grp)
sum_runs  <- function(m) t(rowsum(t(m), group = grp, reorder = FALSE))
mean_runs <- function(m) sweep(sum_runs(m), 2, as.numeric(n_runs[colnames(sum_runs(m))]), "/")
txi2 <- list(counts = sum_runs(txi$counts),
             abundance = mean_runs(txi$abundance),
             length = mean_runs(txi$length),
             countsFromAbundance = "no")

# ---- Sample-level metadata (one row per biological sample) ----
meta <- unique(runs[, c("sample_id", "gsm", "condition",
                        "water_potential_MPa", "replicate")])
rownames(meta) <- meta$sample_id
meta <- meta[colnames(txi2$counts), ]
meta$condition <- factor(make.names(meta$condition))
ref <- "well.watered"
stopifnot(ref %in% levels(meta$condition))
meta$condition <- relevel(meta$condition, ref = ref)
meta$replicate <- factor(meta$replicate)
print(meta)

# ---- DESeq2 ----
dds <- DESeqDataSetFromTximport(txi2, colData = meta, design = ~ condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]   # keep genes with >=10 reads in >=2 samples
dds <- DESeq(dds)
saveRDS(dds, "results/dds.rds")

# ---- Results: each condition vs well-watered ----
for (lv in setdiff(levels(meta$condition), ref)) {
  res <- results(dds, contrast = c("condition", lv, ref), alpha = 0.05)
  cat("\n=== ", lv, " vs ", ref, " ===\n", sep = "")
  summary(res)
  write.csv(as.data.frame(res), paste0("results/deseq2_", lv, "_vs_", ref, ".csv"))
  png(paste0("figures/ma_plot_", lv, ".png"), width = 1800, height = 1500, res = 300)
  plotMA(res, ylim = c(-5, 5), main = paste(lv, "vs", ref))
  dev.off()
}

# ---- Transformed data for QC plots ----
vsd <- vst(dds, blind = TRUE)

# PCA
pca <- plotPCA(vsd, intgroup = c("condition", "replicate"), returnData = TRUE)
pv <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, colour = condition, shape = replicate)) +
  geom_point(size = 3.5) +
  xlab(paste0("PC1: ", pv[1], "% variance")) +
  ylab(paste0("PC2: ", pv[2], "% variance")) +
  theme_bw()
ggsave("figures/pca.png", p, width = 6.5, height = 5, dpi = 300)

# Sample distance heatmap
d <- as.matrix(dist(t(assay(vsd))))
ann <- data.frame(condition = meta$condition, row.names = rownames(meta))
png("figures/sample_distance_heatmap.png", width = 2000, height = 1800, res = 300)
pheatmap(d, annotation_col = ann,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255))
dev.off()
