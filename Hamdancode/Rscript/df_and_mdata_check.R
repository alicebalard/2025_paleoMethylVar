suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

# -----------------------------
# Helpers
# -----------------------------
norm_chr <- function(x) {
  x <- as.character(x)
  ifelse(is.na(x), NA_character_,
         ifelse(str_detect(x, "^chr"), x, paste0("chr", x)))
}

mode1 <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

# ============================================================
# GOAL:
# Find CpGs that are in Atlas "23M" (atlas_annot) AND present in Mdata
# Then (optionally) attach 1-row-per-CpG Atlas annotations to those overlaps
# ============================================================

# -----------------------------
# 1) Atlas keys (23M universe)
#    Make a unique chr:pos key per CpG
# -----------------------------
A_keys <- atlas_annot %>%
  mutate(
    chr_ucsc = norm_chr(chr_ucsc),
    pos_int  = as.integer(pos_int),
    key      = paste0(chr_ucsc, ":", pos_int)
  ) %>%
  filter(!is.na(chr_ucsc), !is.na(pos_int)) %>%
  distinct(key)

# OPTIONAL: collapse atlas_annot to one row per CpG with mode annotation / gene
A_annot_1row <- atlas_annot %>%
  mutate(
    chr_ucsc = norm_chr(chr_ucsc),
    pos_int  = as.integer(pos_int),
    key      = paste0(chr_ucsc, ":", pos_int)
  ) %>%
  filter(!is.na(chr_ucsc), !is.na(pos_int)) %>%
  group_by(key) %>%
  summarise(
    annotation_mode = if ("annotation" %in% names(.)) mode1(annotation) else NA_character_,
    symbol_mode     = if ("SYMBOL" %in% names(.)) mode1(SYMBOL) else NA_character_,
    genename_mode   = if ("GENENAME" %in% names(.)) mode1(GENENAME) else NA_character_,
    mean_distToTSS  = if ("distanceToTSS" %in% names(.)) mean(as.numeric(distanceToTSS), na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

# -----------------------------
# 2) Mdata keys (CpGs present in your sample/data)
#    Tries common schemas: chr/pos, chr_ucsc/pos_int, seqnames/start, cpg "chr:pos", name "chr_pos"
# -----------------------------
M_keys <- Mdata %>%
  mutate(
    chr_tmp = case_when(
      "chr_ucsc" %in% names(.) ~ as.character(chr_ucsc),
      "chr"      %in% names(.) ~ as.character(chr),
      "seqnames" %in% names(.) ~ as.character(seqnames),
      "cpg"      %in% names(.) ~ str_extract(as.character(cpg), "^[^:]+"),
      "name"     %in% names(.) ~ sub("_(.*)", "", as.character(name)),
      TRUE ~ NA_character_
    ),
    pos_tmp = case_when(
      "pos_int" %in% names(.) ~ as.integer(pos_int),
      "pos"     %in% names(.) ~ as.integer(pos),
      "start"   %in% names(.) ~ as.integer(start),
      "cpg"     %in% names(.) ~ as.integer(str_extract(as.character(cpg), "(?<=:)\\d+")),
      "name"    %in% names(.) ~ as.integer(sub(".*_", "", as.character(name))),
      TRUE ~ NA_integer_
    ),
    chr_tmp = norm_chr(chr_tmp),
    pos_tmp = as.integer(pos_tmp),
    key     = paste0(chr_tmp, ":", pos_tmp)
  ) %>%
  filter(!is.na(chr_tmp), !is.na(pos_tmp)) %>%
  distinct(key)

# -----------------------------
# 3) Overlap: CpGs in 23M that are present in Mdata
# -----------------------------
M_in_23M_keys <- inner_join(M_keys, A_keys, by = "key")

cat("Mdata unique CpGs:", nrow(M_keys), "\n")
cat("Atlas 23M unique CpGs:", nrow(A_keys), "\n")
cat("Overlap CpGs (Mdata ∩ 23M):", nrow(M_in_23M_keys), "\n")

# -----------------------------
# 4) OPTIONAL: attach collapsed Atlas annotations to the overlapping CpGs
# -----------------------------
M_in_23M_annot <- M_in_23M_keys %>%
  inner_join(A_annot_1row, by = "key")

cat("Annotated overlap rows:", nrow(M_in_23M_annot), "\n")

# This is your overlap table:
M_in_23M_annot

# OPTIONAL save:
# write.csv(M_in_23M_annot, "Mdata_CpGs_present_in_Atlas23M_with_annotations.csv", row.names = FALSE)
