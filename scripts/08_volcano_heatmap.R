library(DESeq2)
library(ggplot2)
library(pheatmap)

dds <- readRDS("results/dds.rds")
vsd <- vst(dds, blind = TRUE)
ref <- "well.watered"
padj_cut <- 0.05
lfc_cut <- 1

for (lv in c("mild.water.deficit", "severe.water.deficit")) {
  res <- read.csv(paste0("results/deseq2_", lv, "_vs_", ref, ".csv"), row.names = 1)
  res <- res[!is.na(res$padj), ]
  res$status <- "Not significant"
  res$status[res$padj < padj_cut & res$log2FoldChange >= lfc_cut] <- "Up"
  res$status[res$padj < padj_cut & res$log2FoldChange <= -lfc_cut] <- "Down"
  cat("\n", lv, "vs", ref, "\n")
  print(table(res$status))
  write.csv(res[res$status != "Not significant", ],
            paste0("results/sig_", lv, ".csv"))

  p <- ggplot(res, aes(log2FoldChange, -log10(padj), colour = status)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_colour_manual(values = c(Down = "#2166ac",
                                   "Not significant" = "grey70",
                                   Up = "#b2182b")) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
    labs(title = paste(lv, "vs", ref),
         x = "log2 fold change", y = "-log10 adjusted p-value") +
    theme_bw()
  ggsave(paste0("figures/volcano_", lv, ".png"), p, width = 6.5, height = 5, dpi = 300)
}

# Heatmap: top 50 significant genes (severe vs well-watered), by adjusted p-value
sev <- read.csv("results/sig_severe.water.deficit.csv", row.names = 1)
top <- rownames(sev)[order(sev$padj)][1:50]
mat <- assay(vsd)[top, ]
ann <- data.frame(condition = colData(dds)$condition, row.names = colnames(dds))
png("figures/heatmap_top50_severe.png", width = 2000, height = 2400, res = 300)
pheatmap(mat, scale = "row", annotation_col = ann,
         show_rownames = TRUE, fontsize_row = 6)
dev.off()
