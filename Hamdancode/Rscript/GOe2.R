# ============================================================
# FIX: your "top" set is ~half the universe (1594/3082), so ORA has no power.
# Use (A) a smaller TOP gene set for enrichGO + (B) GSEA (gseGO) on ranked genes.
# Self-contained from trail_petrous_annotated.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})

min_cov <- 3
min_app <- 3

# ---- 1) collapse duplicates per (chr,pos,ID), compute beta, apply cov/app thresholds ----
tp_cov <- trail_petrous_annotated %>%
  dplyr::mutate(
    chr = as.character(chr),
    pos = as.integer(pos),
    ID  = as.character(ID),
    cov = as.numeric(cov),
    mthyl = as.numeric(mthyl)
  ) %>%
  dplyr::group_by(chr, pos, ID) %>%
  dplyr::summarise(
    cov   = sum(cov, na.rm = TRUE),
    mthyl = sum(mthyl, na.rm = TRUE),
    beta  = ifelse(cov > 0, mthyl / cov, NA_real_),
    SYMBOL = dplyr::first(na.omit(SYMBOL)),
    .groups = "drop"
  ) %>%
  dplyr::filter(is.finite(beta), cov >= min_cov)

cpg_stats <- tp_cov %>%
  dplyr::group_by(chr, pos) %>%
  dplyr::summarise(
    appearances = n_distinct(ID),
    sd_beta = sd(beta, na.rm = TRUE),
    SYMBOL = paste(sort(unique(na.omit(SYMBOL))), collapse = ";"),
    .groups = "drop"
  ) %>%
  dplyr::filter(appearances >= min_app, is.finite(sd_beta))

# ---- 2) gene-level score: for each gene, take max(sd_beta) across its CpGs ----
gene_scores <- cpg_stats %>%
  dplyr::select(SYMBOL, sd_beta) %>%
  tidyr::separate_rows(SYMBOL, sep = ";") %>%
  dplyr::mutate(SYMBOL = trimws(SYMBOL)) %>%
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(score = max(sd_beta, na.rm = TRUE), .groups = "drop")

# SYMBOL -> ENTREZ
gene_map <- clusterProfiler::bitr(gene_scores$SYMBOL,
                                  fromType = "SYMBOL", toType = "ENTREZID",
                                  OrgDb = org.Hs.eg.db)

gene_scores2 <- gene_scores %>%
  dplyr::inner_join(gene_map, by = "SYMBOL") %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::summarise(score = max(score, na.rm = TRUE), .groups = "drop")

# Universe = all genes that survive thresholds + mapping
bg_entrez <- unique(gene_scores2$ENTREZID)

# Ranked vector for GSEA
geneList <- gene_scores2$score
names(geneList) <- gene_scores2$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

# ---- 3A) ORA with a SMALLER top gene set (this is the key "fix") ----
top_n_genes <- 200  # try 100, 200, 300
top_entrez_small <- names(geneList)[seq_len(min(top_n_genes, length(geneList)))]

ego_bp_small <- enrichGO(
  gene          = top_entrez_small,
  universe      = bg_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.2,
  qvalueCutoff  = 0.2,
  minGSSize     = 3,
  readable      = TRUE
)

# ---- 3B) GSEA (often works even when ORA fails) ----
gse_bp <- gseGO(
  geneList      = geneList,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  minGSSize     = 3,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.2,
  verbose       = FALSE
)

# ---- 4) Plot + save whichever returns results ----
dir.create("GO_fix_outputs", showWarnings = FALSE)

if (nrow(as.data.frame(ego_bp_small)) > 0) {
  write.csv(as.data.frame(ego_bp_small),
            file = "GO_fix_outputs/GO_BP_enrichGO_top200genes_cov3_app3.csv",
            row.names = FALSE)
  
  p_dot <- dotplot(ego_bp_small, showCategory = 20) +
    ggtitle("GO BP (enrichGO): TOP 200 genes by CpG SD (cov>=3, app>=3)")
  ggsave("GO_fix_outputs/GO_BP_dotplot_enrichGO_top200.png", p_dot, width = 10, height = 6, dpi = 300)
  
  p_bar <- barplot(ego_bp_small, showCategory = 20) +
    ggtitle("GO BP (enrichGO): TOP 200 genes by CpG SD (cov>=3, app>=3)")
  ggsave("GO_fix_outputs/GO_BP_barplot_enrichGO_top200.png", p_bar, width = 10, height = 6, dpi = 300)
  
  print(p_dot)
} else {
  message("enrichGO still returned 0 terms with TOP ", top_n_genes, " genes. (Try top_n_genes = 100 or 300.)")
}

if (nrow(as.data.frame(gse_bp)) > 0) {
  write.csv(as.data.frame(gse_bp),
            file = "GO_fix_outputs/GO_BP_gseGO_ranked_cov3_app3.csv",
            row.names = FALSE)
  
  p_gse <- dotplot(gse_bp, showCategory = 20) +
    ggtitle("GO BP (gseGO): Ranked genes by CpG SD (cov>=3, app>=3)")
  ggsave("GO_fix_outputs/GO_BP_dotplot_gseGO.png", p_gse, width = 10, height = 6, dpi = 300)
  
  print(p_gse)
} else {
  message("gseGO returned 0 terms at p/q<=0.2. Try pvalueCutoff=0.5 or check score distribution.")
}

message("Done. Check the folder: GO_fix_outputs/")













# ============================================================
# PETROUS: Robust GO GSEA (gseGO) + ORA fallback
# Fixes the "missing value where TRUE/FALSE needed" error by:
#   - building a clean, named, finite, sorted geneList (ENTREZID -> score)
#   - deduplicating ENTREZIDs
#   - breaking ties slightly
# Outputs:
#   - petrous_GSEA_GO_BP_lenient.csv
#   - petrous_GSEA_GO_BP_strict.csv
#   - petrous_ORA_GO_BP_lenient.csv
#   - plots/petrous_GSEA_GO_BP_dotplot.png (if results)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ---- Bioconductor installs (only if missing) ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, update = FALSE, ask = FALSE)
  }
}
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
})

dir.create("plots", showWarnings = FALSE)

# -----------------------
# REQUIREMENTS
# -----------------------
# You must already have:
#   petrous_var   : per-CpG summary table containing SYMBOL_mode and a variability score
#                  (columns: SYMBOL_mode + either sd_beta OR var_score)
#   top_var_cpgs  : (optional) for ORA foreground; if not present we create from petrous_var
#
# If you only have petrous_var, that's enough for GSEA and ORA.

# -----------------------
# Choose score column
# -----------------------
score_col <- if ("var_score" %in% names(petrous_var)) "var_score" else "sd_beta"
if (!score_col %in% names(petrous_var)) stop("petrous_var must contain 'sd_beta' or 'var_score'.")

# -----------------------
# 1) Build a clean geneList for GSEA
#    Rank genes by max CpG variability assigned to that gene
# -----------------------
gene_rank <- petrous_var %>%
  dplyr::filter(!is.na(SYMBOL_mode), SYMBOL_mode != "") %>%
  dplyr::group_by(SYMBOL_mode) %>%
  dplyr::summarise(score = max(.data[[score_col]], na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(is.finite(score))

cat("Unique SYMBOLs with finite score:", nrow(gene_rank), "\n")

map <- suppressMessages(
  bitr(gene_rank$SYMBOL_mode, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
)

gene_rank2 <- gene_rank %>%
  dplyr::inner_join(map, by = c("SYMBOL_mode" = "SYMBOL")) %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::summarise(score = max(score, na.rm = TRUE), .groups = "drop") %>%  # dedupe ENTREZ
  dplyr::filter(is.finite(score))

cat("Unique ENTREZIDs with finite score:", nrow(gene_rank2), "\n")
if (nrow(gene_rank2) < 20) warning("Very few mapped ENTREZIDs; GO results may be empty/unstable.")

geneList <- gene_rank2$score
names(geneList) <- gene_rank2$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

# break exact ties slightly (helps stability when many ties)
set.seed(1)
geneList <- geneList + runif(length(geneList), -1e-10, 1e-10)

# final sanity checks
stopifnot(is.numeric(geneList), !anyNA(geneList), all(is.finite(geneList)), any(names(geneList) != ""))
cat("geneList ready. Range:", paste(range(geneList), collapse=" .. "), "\n")

# -----------------------
# 2) GSEA (lenient then stricter)
# -----------------------
gse_bp_lenient <- gseGO(
  geneList      = geneList,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  minGSSize     = 10,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,
  verbose       = FALSE,
  eps           = 0
)

gse_df_lenient <- as.data.frame(gse_bp_lenient)
cat("GSEA BP terms (lenient):", nrow(gse_df_lenient), "\n")
write.csv(gse_df_lenient, "petrous_GSEA_GO_BP_lenient.csv", row.names = FALSE)

gse_bp_strict <- gseGO(
  geneList      = geneList,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  minGSSize     = 10,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.2,
  verbose       = FALSE,
  eps           = 0
)

gse_df_strict <- as.data.frame(gse_bp_strict)
cat("GSEA BP terms (strict):", nrow(gse_df_strict), "\n")
write.csv(gse_df_strict, "petrous_GSEA_GO_BP_strict.csv", row.names = FALSE)

if (nrow(gse_df_lenient) > 0) {
  p_gsea <- dotplot(gse_bp_lenient, showCategory = 20) +
    ggtitle("Petrous GO BP GSEA (lenient)")
  ggsave("plots/petrous_GSEA_GO_BP_dotplot.png", p_gsea, width = 9, height = 6, dpi = 300)
}

# -----------------------
# 3) ORA fallback (top genes vs background genes)
#    Foreground: genes from the top variable CpGs
#    Background: genes from all filtered CpGs (petrous_var)
# -----------------------
# Define top_var_cpgs if not already created
if (!exists("top_var_cpgs")) {
  top_var_cpgs <- petrous_var %>% dplyr::arrange(dplyr::desc(.data[[score_col]])) %>% dplyr::slice_head(n = 200)
}

genes_top <- top_var_cpgs %>%
  dplyr::transmute(SYMBOL = SYMBOL_mode) %>%
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
  dplyr::distinct() %>%
  dplyr::pull(SYMBOL)

genes_bg <- petrous_var %>%
  dplyr::transmute(SYMBOL = SYMBOL_mode) %>%
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
  dplyr::distinct() %>%
  dplyr::pull(SYMBOL)

top_map2 <- suppressMessages(bitr(genes_top, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db))
bg_map2  <- suppressMessages(bitr(genes_bg,  fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db))

top_entrez <- unique(top_map2$ENTREZID)
bg_entrez  <- unique(bg_map2$ENTREZID)

ego_bp_lenient <- enrichGO(
  gene          = top_entrez,
  universe      = bg_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 1,
  qvalueCutoff  = 1,
  minGSSize     = 5,
  readable      = TRUE
)

ego_df_lenient <- as.data.frame(ego_bp_lenient)
cat("ORA BP terms (lenient):", nrow(ego_df_lenient), "\n")
write.csv(ego_df_lenient, "petrous_ORA_GO_BP_lenient.csv", row.names = FALSE)

cat("Done. Files written:\n",
    " - petrous_GSEA_GO_BP_lenient.csv\n",
    " - petrous_GSEA_GO_BP_strict.csv\n",
    " - petrous_ORA_GO_BP_lenient.csv\n",
    " - plots/petrous_GSEA_GO_BP_dotplot.png (if any GSEA terms)\n")












