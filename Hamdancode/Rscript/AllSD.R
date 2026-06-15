# ============================================================
# FAST overlap-only pipeline (subset trail to combinedSD CpGs first)
# Same structure as before: prep -> subset -> collapse -> filter -> SD -> overlap -> corr -> plot
# ============================================================

combinedSD <- read.csv("newset/mean_cpg_std_across_all_datasets.csv")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

min_cov    <- 3
min_appear <- 3

# ---------------------------
# 0) combinedSD prep  (replaces neuro)
# ---------------------------
combo <- as.data.table(combinedSD)

combo[, CpG_Name := as.character(CpG_Name)]
combo[, CpG_Name := gsub(":", "_", CpG_Name)]  # just in case any are chr:pos
combo[, SD_combo := as.numeric(Mean_Standard_Deviation)]
combo <- combo[is.finite(SD_combo), .(CpG_Name, SD_combo)]

# make chrpos to match trail (chr:pos)
combo[, chrpos := gsub("_", ":", CpG_Name)]

# ---------------------------
# 1) trail prep
# ---------------------------
trail <- as.data.table(trail_petrous2)
trail[, `:=`(
  chr    = as.character(chr),
  pos    = as.integer(pos),
  ID     = as.character(ID),
  cov    = as.numeric(cov),
  mthyl  = as.numeric(mthyl),
  chrpos = as.character(chrpos)
)]

# 1) SUBSET to combinedSD CpGs (massive speedup)
trail_sub <- trail[chrpos %chin% combo$chrpos]

# 2) Collapse duplicates to 1 row per (chrpos, ID)
trail_per_id <- trail_sub[, .(
  cov   = sum(cov,  na.rm = TRUE),
  mthyl = sum(mthyl, na.rm = TRUE)
), by = .(chrpos, chr, pos, ID)]

trail_per_id[, beta := fifelse(cov > 0, mthyl / cov, NA_real_)]

# 3) Filters: cov>=3 (ID-level), appearances>=3 (IDs per CpG)
trail_filt <- trail_per_id[
  cov >= min_cov & is.finite(beta)
][
  , if (.N >= min_appear) .SD, by = .(chrpos, chr, pos)
]

# 4) SD(beta) per CpG
trail_sd <- trail_filt[, .(
  n_ids   = .N,                       # per_id is unique per ID now
  sd_beta = sd(beta, na.rm = TRUE)
), by = .(chrpos, chr, pos)]

trail_sd <- trail_sd[is.finite(sd_beta)]
trail_sd[, CpG_Name := gsub(":", "_", chrpos)]

# 5) Overlap + correlation (trail_sd vs combinedSD)
overlap <- merge(
  trail_sd[, .(CpG_Name, sd_beta, n_ids)],
  combo[, .(CpG_Name, SD_combo)],
  by = "CpG_Name"
)

cat("Trail rows (subset):", nrow(trail_sub), "\n")
cat("CpGs with SD(beta):", nrow(trail_sd), "\n")
cat("Overlapping CpGs:", nrow(overlap), "\n")

pearson_res  <- cor.test(overlap$sd_beta, overlap$SD_combo, method = "pearson")
spearman_res <- cor.test(overlap$sd_beta, overlap$SD_combo, method = "spearman")

r_txt   <- sprintf("Pearson r = %.3f (p=%.2g)",  unname(pearson_res$estimate),  pearson_res$p.value)
rho_txt <- sprintf("Spearman ρ = %.3f (p=%.2g)", unname(spearman_res$estimate), spearman_res$p.value)

cat("\n", r_txt, "\n", rho_txt, "\n", sep="")

                            
# 6) Plot
x_annot <- as.numeric(quantile(overlap$sd_beta, 0.05, na.rm = TRUE))
y_annot <- as.numeric(quantile(overlap$SD_combo, 0.95, na.rm = TRUE))

p <- ggplot(overlap, aes(x = sd_beta, y = SD_combo)) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  annotate("text", x = x_annot, y = y_annot, hjust = 0, vjust = -5,
           label = paste(r_txt, rho_txt, sep = "\n")) +
  labs(
    title = paste0("Overlap CpGs (cov≥", min_cov, ", IDs≥", min_appear, "): Trail vs CombinedSD"),
    x = "Trail SD(beta) across IDs",
    y = "Combined SD"
  ) +
  theme_bw()

print(p)
ggsave("trail_combinedSD_scatter.png", p, width = 7, height = 5, dpi = 300)                         





















# ============================================================
# Heatmap (like your screenshot) for: Trail SD(beta) vs CombinedSD (Mean SD)
# Grid over min_cov (ID-level) x min_appear (IDs per CpG)
# ============================================================

combinedSD <- read.csv("newset/mean_cpg_std_across_all_datasets.csv")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---- grid (same as before / like the example) ----
cov_grid <- c(1, 2, 3, 5, 10)
app_grid <- c(2, 3, 5, 8, 10)
cor_method <- "spearman"

# ---------------------------
# 0) combinedSD prep
# ---------------------------
combo_dt <- as.data.table(combinedSD)
combo_dt[, CpG_Name := gsub(":", "_", as.character(CpG_Name))]
combo_dt[, SD_combo := as.numeric(Mean_Standard_Deviation)]
combo_dt <- combo_dt[is.finite(SD_combo), .(CpG_Name, SD_combo)]
combo_dt[, chrpos := gsub("_", ":", CpG_Name)]
setkey(combo_dt, chrpos)

# ---------------------------
# 1) trail prep
# ---------------------------
trail <- as.data.table(trail_petrous2)
trail[, `:=`(
  chr  = as.character(chr),
  pos  = as.integer(pos),
  ID   = as.character(ID),
  cov  = as.numeric(cov),
  mthyl= as.numeric(mthyl)
)]

# make chrpos if missing
if (!("chrpos" %in% names(trail))) {
  trail[, chrpos := paste0(chr, ":", pos)]
}
trail[, chrpos := as.character(chrpos)]

# 1) SUBSET to combinedSD CpGs (massive speedup)
trail_sub <- trail[chrpos %chin% combo_dt$chrpos]

# 2) Collapse duplicates -> 1 row per (chrpos, ID)
trail_per_id <- trail_sub[, .(
  cov   = sum(cov,  na.rm = TRUE),
  mthyl = sum(mthyl,na.rm = TRUE)
), by = .(chrpos, ID)]

trail_per_id[, beta := fifelse(cov > 0, mthyl / cov, NA_real_)]
trail_per_id <- trail_per_id[is.finite(beta)]

# ---------------------------
# helper: safe cor
# ---------------------------
safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3) return(list(r = NA_real_, p = NA_real_))
  out <- tryCatch({
    ct <- cor.test(x, y, method = method)
    list(r = as.numeric(unname(ct$estimate)), p = as.numeric(ct$p.value))
  }, error = function(e) list(r = NA_real_, p = NA_real_))
  out
}

# ---------------------------
# 2) Compute heatmap table
# ---------------------------
res_list <- vector("list", length(cov_grid) * length(app_grid))
k <- 1L

for (min_cov in cov_grid) {
  
  dt_cov <- trail_per_id[cov >= min_cov]
  
  # SD(beta) + n IDs per CpG
  stats_cov <- dt_cov[, .(
    n_ids   = .N,
    sd_beta = sd(beta, na.rm = TRUE)
  ), by = .(chrpos)]
  
  stats_cov <- stats_cov[is.finite(sd_beta)]
  
  for (min_app in app_grid) {
    
    keep <- stats_cov[n_ids >= min_app]
    
    ov <- merge(
      keep[, .(chrpos, sd_beta)],
      combo_dt[, .(chrpos, SD_combo)],
      by = "chrpos"
    )
    
    n_ov <- nrow(ov)
    ct <- safe_cor(ov$sd_beta, ov$SD_combo, method = cor_method)
    
    res_list[[k]] <- data.table(
      min_cov    = min_cov,
      min_appear = min_app,
      cor        = ct$r,
      p_value    = ct$p,
      n_cpg      = n_ov
    )
    k <- k + 1L
  }
}

res <- rbindlist(res_list)

# order axes
res[, min_cov := factor(min_cov, levels = cov_grid)]
res[, min_appear := factor(min_appear, levels = app_grid)]

# labels in tiles (cor only, like your screenshot)
res[, lab := ifelse(is.na(cor), "NA", sprintf("%.2f", cor))]
# if you ALSO want counts per square:
# res[, lab := ifelse(is.na(cor), "NA", sprintf("%.2f\n(n=%d)", cor, n_cpg))]

# ---------------------------
# 3) Plot heatmap
# ---------------------------
p_heat <- ggplot(res, aes(x = min_appear, y = min_cov, fill = cor)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = lab), size = 4) +
  scale_fill_gradient2(
    name = "cor",
    midpoint = 0,
    low = "#0B1D39", mid = "#2C6AA3", high = "#5DB8FF",
    na.value = "grey90"
  ) +
  labs(
    title = paste0("Method: ", cor_method),
    x = "number of inds (>=X)",
    y = "Minimum coverage (per ID at CpG)"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())

print(p_heat)
ggsave("trail_combinedSD_heatmap.png", p_heat, width = 8, height = 6, dpi = 300)

# optional: save the table behind the heatmap
fwrite(res, "trail_combinedSD_heatmap_table.csv")












suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---- grid ----
cov_grid <- c(1, 2, 3, 5, 10)
app_grid <- c(2, 3, 5, 8, 10)
cor_method <- "pearson"   # <-- changed

# ---------------------------
# 0) combinedSD prep
# ---------------------------
combo_dt <- as.data.table(combinedSD)
combo_dt[, CpG_Name := gsub(":", "_", as.character(CpG_Name))]
combo_dt[, SD_combo := as.numeric(Mean_Standard_Deviation)]
combo_dt <- combo_dt[is.finite(SD_combo), .(CpG_Name, SD_combo)]
combo_dt[, chrpos := gsub("_", ":", CpG_Name)]
setkey(combo_dt, chrpos)

# ---------------------------
# 1) trail prep
# ---------------------------
trail <- as.data.table(trail_petrous2)
trail[, `:=`(
  chr  = as.character(chr),
  pos  = as.integer(pos),
  ID   = as.character(ID),
  cov  = as.numeric(cov),
  mthyl= as.numeric(mthyl)
)]

# make chrpos if missing
if (!("chrpos" %in% names(trail))) {
  trail[, chrpos := paste0(chr, ":", pos)]
}
trail[, chrpos := as.character(chrpos)]

# subset to combinedSD CpGs (speedup)
trail_sub <- trail[chrpos %chin% combo_dt$chrpos]

# collapse duplicates -> 1 row per (chrpos, ID)
trail_per_id <- trail_sub[, .(
  cov   = sum(cov,  na.rm = TRUE),
  mthyl = sum(mthyl,na.rm = TRUE)
), by = .(chrpos, ID)]

trail_per_id[, beta := fifelse(cov > 0, mthyl / cov, NA_real_)]
trail_per_id <- trail_per_id[is.finite(beta)]

# ---------------------------
# helper: safe cor
# ---------------------------
safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3) return(list(r = NA_real_, p = NA_real_))
  out <- tryCatch({
    ct <- cor.test(x, y, method = method)
    list(r = as.numeric(unname(ct$estimate)), p = as.numeric(ct$p.value))
  }, error = function(e) list(r = NA_real_, p = NA_real_))
  out
}

# ---------------------------
# 2) Compute heatmap table
# ---------------------------
res_list <- vector("list", length(cov_grid) * length(app_grid))
k <- 1L

for (min_cov in cov_grid) {
  
  dt_cov <- trail_per_id[cov >= min_cov]
  
  stats_cov <- dt_cov[, .(
    n_ids   = .N,
    sd_beta = sd(beta, na.rm = TRUE)
  ), by = .(chrpos)]
  
  stats_cov <- stats_cov[is.finite(sd_beta)]
  
  for (min_app in app_grid) {
    
    keep <- stats_cov[n_ids >= min_app]
    
    ov <- merge(
      keep[, .(chrpos, sd_beta)],
      combo_dt[, .(chrpos, SD_combo)],
      by = "chrpos"
    )
    
    n_ov <- nrow(ov)
    ct <- safe_cor(ov$sd_beta, ov$SD_combo, method = cor_method)
    
    res_list[[k]] <- data.table(
      min_cov    = min_cov,
      min_appear = min_app,
      cor        = ct$r,
      p_value    = ct$p,
      n_cpg      = n_ov
    )
    k <- k + 1L
  }
}

res <- rbindlist(res_list)

# order axes
res[, min_cov := factor(min_cov, levels = cov_grid)]
res[, min_appear := factor(min_appear, levels = app_grid)]

# label: correlation + n CpGs in each tile
res[, lab := ifelse(is.na(cor), sprintf("NA\n(n=%d)", n_cpg),
                    sprintf("%.2f\n(n=%d)", cor, n_cpg))]

# ---------------------------
# 3) Plot heatmap
# ---------------------------
p_heat <- ggplot(res, aes(x = min_appear, y = min_cov, fill = cor)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = lab), size = 3.6, lineheight = 0.95) +
  scale_fill_gradient2(
    name = "cor",
    midpoint = 0,
    low = "#0B1D39", mid = "#2C6AA3", high = "#5DB8FF",
    na.value = "grey90"
  ) +
  labs(
    title = paste0("Atlas and ancient standard deviation correlation at different thresholds. Method: ", cor_method),
    x = "Minimum appearances (distinct IDs per CpG)",
    y = "Minimum coverage (per ID at CpG)"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())

print(p_heat)
ggsave("trail_combinedSD_heatmap_pearson.png", p_heat, width = 8, height = 6, dpi = 300)

# optional: save the table behind the heatmap
fwrite(res, "trail_combinedSD_heatmap_table_pearson.csv")









suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---- grid ----
cov_grid <- c(1, 2, 3, 5, 10)
app_grid <- c(2, 3, 5, 8, 10)
cor_method <- "spearman"

# ---------------------------
# 0) combinedSD prep
# ---------------------------
combo_dt <- as.data.table(combinedSD)
combo_dt[, CpG_Name := gsub(":", "_", as.character(CpG_Name))]
combo_dt[, SD_combo := as.numeric(Mean_Standard_Deviation)]
combo_dt <- combo_dt[is.finite(SD_combo), .(CpG_Name, SD_combo)]
combo_dt[, chrpos := gsub("_", ":", CpG_Name)]
setkey(combo_dt, chrpos)

# ---------------------------
# 1) trail prep
# ---------------------------
trail <- as.data.table(trail_petrous2)
trail[, `:=`(
  chr   = as.character(chr),
  pos   = as.integer(pos),
  ID    = as.character(ID),
  cov   = as.numeric(cov),
  mthyl = as.numeric(mthyl)
)]

if (!("chrpos" %in% names(trail))) {
  trail[, chrpos := paste0(chr, ":", pos)]
}
trail[, chrpos := as.character(chrpos)]

trail_sub <- trail[chrpos %chin% combo_dt$chrpos]

trail_per_id <- trail_sub[, .(
  cov   = sum(cov, na.rm = TRUE),
  mthyl = sum(mthyl, na.rm = TRUE)
), by = .(chrpos, ID)]

trail_per_id[, beta := fifelse(cov > 0, mthyl / cov, NA_real_)]
trail_per_id <- trail_per_id[is.finite(beta)]

# ---------------------------
# helper: safe cor
# ---------------------------
safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3) return(list(r = NA_real_, p = NA_real_))
  
  tryCatch({
    ct <- cor.test(x, y, method = method)
    list(r = as.numeric(unname(ct$estimate)), p = as.numeric(ct$p.value))
  }, error = function(e) {
    list(r = NA_real_, p = NA_real_)
  })
}

# ---------------------------
# 2) Compute heatmap table
# ---------------------------
res_list <- vector("list", length(cov_grid) * length(app_grid))
k <- 1L

for (min_cov in cov_grid) {
  dt_cov <- trail_per_id[cov >= min_cov]
  
  stats_cov <- dt_cov[, .(
    n_ids   = .N,
    sd_beta = sd(beta, na.rm = TRUE)
  ), by = .(chrpos)]
  
  stats_cov <- stats_cov[is.finite(sd_beta)]
  
  for (min_app in app_grid) {
    keep <- stats_cov[n_ids >= min_app]
    
    ov <- merge(
      keep[, .(chrpos, sd_beta)],
      combo_dt[, .(chrpos, SD_combo)],
      by = "chrpos"
    )
    
    n_ov <- nrow(ov)
    ct <- safe_cor(ov$sd_beta, ov$SD_combo, method = cor_method)
    
    res_list[[k]] <- data.table(
      min_cov    = min_cov,
      min_appear = min_app,
      cor        = ct$r,
      p_value    = ct$p,
      n_cpg      = n_ov
    )
    k <- k + 1L
  }
}

res <- rbindlist(res_list)

res[, min_cov := factor(min_cov, levels = cov_grid)]
res[, min_appear := factor(min_appear, levels = app_grid)]

# correlation + number of CpGs in each tile
res[, lab := ifelse(
  is.na(cor),
  sprintf("NA\n(n=%d)", n_cpg),
  sprintf("%.2f\n(n=%d)", cor, n_cpg)
)]

# ---------------------------
# 3) Plot heatmap
# ---------------------------
p_heat <- ggplot(res, aes(x = min_appear, y = min_cov, fill = cor)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = lab), size = 3.6, lineheight = 0.95) +
  scale_fill_gradient2(
    name = "rho",
    midpoint = 0,
    low = "#0B1D39", mid = "#2C6AA3", high = "#5DB8FF",
    na.value = "grey90"
  ) +
  labs(
    title = paste0(
      "Atlas and ancient standard deviation correlation at different thresholds. Method: ",
      cor_method
    ),
    x = "Minimum appearances (distinct IDs per CpG)",
    y = "Minimum coverage (per ID at CpG)"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())

print(p_heat)
ggsave("trail_combinedSD_heatmap_spearman.png", p_heat, width = 8, height = 6, dpi = 300)
fwrite(res, "trail_combinedSD_heatmap_table_spearman.csv")
