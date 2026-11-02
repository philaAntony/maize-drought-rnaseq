# Gene lists (up/down) and background for GO enrichment
dir.create("enrichment", showWarnings = FALSE)
for (lv in c("mild.water.deficit", "severe.water.deficit")) {
  s <- read.csv(paste0("results/sig_", lv, ".csv"), row.names = 1)
  writeLines(rownames(s)[s$log2FoldChange > 0], paste0("enrichment/", lv, "_up.txt"))
  writeLines(rownames(s)[s$log2FoldChange < 0], paste0("enrichment/", lv, "_down.txt"))
}
r <- read.csv("results/deseq2_severe.water.deficit_vs_well.watered.csv", row.names = 1)
writeLines(rownames(r)[!is.na(r$padj)], "enrichment/universe.txt")
