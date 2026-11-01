suppressMessages({
  library(clusterProfiler); library(GO.db)
  library(AnnotationDbi); library(ggplot2)
})

# Source: MaizeGDB Pannzer GO predictions for B73 v5, downloaded with
# wget https://download.maizegdb.org/GeneFunction_and_Expression/Pannzer_GO_Terms/B73.out \
#   -O reference/pannzer_B73.out
ppv_cut <- 0.5

# ---- 1. Gene-to-GO table ----
pz <- read.delim("reference/pannzer_B73.out", sep = "\t", header = TRUE, quote = "",
                 comment.char = "", colClasses = "character", fill = TRUE)
go <- pz[pz$type %in% c("BP_ARGOT", "MF_ARGOT", "CC_ARGOT"), ]
go$PPV <- as.numeric(go$PPV)
go <- go[!is.na(go$PPV) & go$PPV >= ppv_cut, ]
go$gene <- sub("_P[0-9]+$", "", go$qpid)
go$term <- paste0("GO:", go$id)
gene2go <- unique(go[, c("gene", "term")])

# ---- 2. Add parent terms (true-path rule) ----
anc <- c(as.list(GOBPANCESTOR), as.list(GOMFANCESTOR), as.list(GOCCANCESTOR))
tt <- unique(gene2go$term)
exp <- lapply(tt, function(t) c(t, anc[[t]]))
map <- data.frame(term = rep(tt, lengths(exp)), anc = unlist(exp, use.names = FALSE),
                  stringsAsFactors = FALSE)
full <- unique(merge(gene2go, map, by = "term")[, c("anc", "gene")])
colnames(full) <- c("term", "gene")
full <- full[full$term != "all", ]

info <- AnnotationDbi::select(GO.db, keys = unique(full$term),
                              columns = c("TERM", "ONTOLOGY"), keytype = "GOID")
info <- info[!is.na(info$TERM), ]
full <- full[full$term %in% info$GOID, ]

# ---- 3. Background ----
univ <- readLines("enrichment/universe.txt")
cat("Tested genes with at least one GO term:",
    length(intersect(univ, unique(full$gene))), "of", length(univ), "\n")

# ---- 4. Enrichment ----
lists <- c("mild.water.deficit_up", "mild.water.deficit_down",
           "severe.water.deficit_up", "severe.water.deficit_down")
all_res <- list()
for (l in lists) {
  g <- readLines(paste0("enrichment/", l, ".txt"))
  for (ont in c("BP", "MF", "CC")) {
    t2n <- info[info$ONTOLOGY == ont, c("GOID", "TERM")]
    t2g <- full[full$term %in% t2n$GOID, c("term", "gene")]
    e <- enricher(g, universe = univ, TERM2GENE = t2g, TERM2NAME = t2n,
                  pvalueCutoff = 1, qvalueCutoff = 1,
                  minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH")
    if (!is.null(e) && nrow(e@result) > 0) {
      r <- e@result
      r$list <- l
      r$ontology <- ont
      all_res[[paste(l, ont)]] <- r
    }
  }
}

res <- do.call(rbind, all_res)
sig <- res[res$p.adjust < 0.05, ]
sig$GeneRatioNum <- sapply(strsplit(sig$GeneRatio, "/"),
                           function(x) as.numeric(x[1]) / as.numeric(x[2]))
write.csv(sig, "results/go_enrichment_significant.csv", row.names = FALSE)
cat("\nSignificant terms (adjusted p < 0.05):\n")
print(table(sig$list, sig$ontology))

# ---- 5. Dot plots (biological process) ----
for (l in lists) {
  d <- sig[sig$list == l & sig$ontology == "BP", ]
  if (nrow(d) == 0) next
  d <- head(d[order(d$p.adjust), ], 15)
  d$Description <- factor(d$Description, levels = rev(unique(d$Description)))
  p <- ggplot(d, aes(GeneRatioNum, Description, size = Count, colour = p.adjust)) +
    geom_point() +
    labs(title = paste("GO biological process:", l),
         x = "Gene ratio", y = NULL) +
    theme_bw()
  ggsave(paste0("figures/go_bp_", l, ".png"), p, width = 8, height = 5.5, dpi = 300)
}
