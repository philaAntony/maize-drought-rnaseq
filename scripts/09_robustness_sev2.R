suppressMessages({library(DESeq2); library(ggplot2)})

dds_full <- readRDS("results/dds.rds")
ref <- "well.watered"
lv <- "severe.water.deficit"
padj_cut <- 0.05
lfc_cut <- 1
cd <- colData(dds_full)[, c("sample_id", "condition", "replicate")]

# Refit the model with some samples dropped (same pipeline as 07_deseq2.R)
fit <- function(drop = character(0)) {
  keep <- setdiff(colnames(dds_full), drop)
  d <- DESeqDataSetFromMatrix(countData = round(counts(dds_full)[, keep]),
                              colData = cd[keep, ], design = ~ condition)
  assays(d)[["avgTxLength"]] <- assays(dds_full)[["avgTxLength"]][, keep]
  d <- d[rowSums(counts(d) >= 10) >= 2, ]
  d <- DESeq(d, quiet = TRUE)
  as.data.frame(results(d, contrast = c("condition", lv, ref), alpha = 0.05))
}

# Genes passing padj < 0.05 and |log2FC| >= 1
sig <- function(res) {
  ok <- !is.na(res$padj)
  list(up = rownames(res)[ok & res$padj < padj_cut & res$log2FoldChange >= lfc_cut],
       down = rownames(res)[ok & res$padj < padj_cut & res$log2FoldChange <= -lfc_cut])
}

full <- fit()
s_full <- sig(full)
cat("All 9 samples: up", length(s_full$up), "down", length(s_full$down), "\n")

drops <- c("SEV_1", "SEV_2", "SEV_3")
fits <- setNames(lapply(drops, fit), drops)

out <- do.call(rbind, lapply(drops, function(x) {
  r <- fits[[x]]
  s <- sig(r)
  common <- intersect(rownames(full), rownames(r))
  data.frame(
    dropped = x,
    up = length(s$up),
    down = length(s$down),
    pct_up_kept = round(100 * mean(s_full$up %in% s$up), 1),
    pct_down_kept = round(100 * mean(s_full$down %in% s$down), 1),
    new_up = sum(!s$up %in% s_full$up),
    new_down = sum(!s$down %in% s_full$down),
    lfc_spearman = round(cor(full[common, "log2FoldChange"],
                             r[common, "log2FoldChange"],
                             method = "spearman", use = "complete.obs"), 3)
  )
}))
print(out)
write.csv(out, "results/robustness_leave_one_out.csv", row.names = FALSE)

# Fold changes with and without SEV_2
r2 <- fits[["SEV_2"]]
common <- intersect(rownames(full), rownames(r2))
df <- data.frame(all9 = full[common, "log2FoldChange"],
                 no_sev2 = r2[common, "log2FoldChange"])
p <- ggplot(df, aes(all9, no_sev2)) +
  geom_point(alpha = 0.3, size = 0.6) +
  geom_abline(slope = 1, intercept = 0, colour = "red") +
  labs(x = "log2 fold change (all 9 samples)",
       y = "log2 fold change (SEV_2 removed)",
       title = "Severe vs well-watered") +
  theme_bw()
ggsave("figures/robustness_sev2_lfc.png", p, width = 5, height = 5, dpi = 300)
