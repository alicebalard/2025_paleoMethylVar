neuro <- read.csv("newset/cpg_std_values_Brain-Neuronal.csv")
getwd()
setwd("/Users/hamda/Desktop/UGI_thingy/3rdyear/GhubAlice/2025_paleoMethylVar/Hamdancode")


head(neuro)

trail_petrous2

# ============================================================
# FAST overlap-only pipeline (subset trail to neuro CpGs first)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

min_cov    <- 3
min_appear <- 3

# ---- neuro prep ----
neuro <- as.data.table(neuro)
neuro[, CpG_Name := as.character(CpG_Name)]
neuro[, CpG_Name := gsub(":", "_", CpG_Name)]
neuro[, SD_neuro := as.numeric(Standard_Deviation)]
neuro <- neuro[is.finite(SD_neuro), .(CpG_Name, SD_neuro)]

# make chrpos to match trail (chr:pos)
neuro[, chrpos := gsub("_", ":", CpG_Name)]

# ---- trail prep ----
trail <- as.data.table(trail_petrous2)
trail[, `:=`(
  chr   = as.character(chr),
  pos   = as.integer(pos),
  ID    = as.character(ID),
  cov   = as.numeric(cov),
  mthyl = as.numeric(mthyl),
  chrpos = as.character(chrpos)
)]

# 1) SUBSET to neuro CpGs (massive speedup)
trail_sub <- trail[chrpos %chin% neuro$chrpos]

# 2) Collapse duplicates to 1 row per (chrpos, ID)
trail_per_id <- trail_sub[, .(
  cov  = sum(cov,  na.rm = TRUE),
  mthyl= sum(mthyl,na.rm = TRUE)
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
  n_ids   = .N,                       # because per_id is unique per ID now
  sd_beta = sd(beta, na.rm = TRUE)
), by = .(chrpos, chr, pos)]

trail_sd <- trail_sd[is.finite(sd_beta)]
trail_sd[, CpG_Name := gsub(":", "_", chrpos)]

# 5) Overlap + correlation
overlap <- merge(trail_sd[, .(CpG_Name, sd_beta, n_ids)],
                 neuro[, .(CpG_Name, SD_neuro)],
                 by = "CpG_Name")

cat("Trail rows (subset):", nrow(trail_sub), "\n")
cat("CpGs with SD(beta):", nrow(trail_sd), "\n")
cat("Overlapping CpGs:", nrow(overlap), "\n")

pearson_res  <- cor.test(overlap$sd_beta, overlap$SD_neuro, method = "pearson")
spearman_res <- cor.test(overlap$sd_beta, overlap$SD_neuro, method = "spearman")

r_txt   <- sprintf("Pearson r = %.3f (p=%.2g)",  unname(pearson_res$estimate),  pearson_res$p.value)
rho_txt <- sprintf("Spearman ρ = %.3f (p=%.2g)", unname(spearman_res$estimate), spearman_res$p.value)

cat("\n", r_txt, "\n", rho_txt, "\n", sep="")

# 6) Plot
x_annot <- as.numeric(quantile(overlap$sd_beta, 0.05, na.rm = TRUE))
y_annot <- as.numeric(quantile(overlap$SD_neuro, 0.95, na.rm = TRUE))

p <- ggplot(overlap, aes(x = sd_beta, y = SD_neuro)) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Ancient petrous SD",
    y = "Modern Atlas neuronal SD",
    title = NULL,
    subtitle = NULL
  ) +
  coord_cartesian(xlim = c(0, 0.375), ylim = c(0, 0.10)) +
  theme_bw()

print(p)

cat("n =", n_overlap, "\n")
cat("Spearman r =", round(spearman_r, 3), "\n")

ggsave("trail_neuro_sd_scatter.png", p, width = 7, height = 5, dpi = 300)




 

# ============================================================
# Heatmap (like your screenshot) of Spearman correlation across:
#   min_cov (ID-level coverage)  x  min_appear (IDs per CpG)
# Comparing: Trail SD(beta) vs Neuro SD
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---- grid (matches your example image) ----
cov_grid <- c(1, 2, 3, 5, 10)
app_grid <- c(2, 3, 5, 8, 10)
cor_method <- "spearman"

# ---- neuro prep ----
neuro_dt <- as.data.table(neuro)
neuro_dt[, CpG_Name := gsub(":", "_", as.character(CpG_Name))]
neuro_dt[, SD_neuro := as.numeric(Standard_Deviation)]
neuro_dt <- neuro_dt[is.finite(SD_neuro), .(CpG_Name, SD_neuro)]
neuro_dt[, chrpos := gsub("_", ":", CpG_Name)]
setkey(neuro_dt, chrpos)

# ---- trail prep ----
trail <- as.data.table(trail_petrous2)
trail[, `:=`(
  chr   = as.character(chr),
  pos   = as.integer(pos),
  ID    = as.character(ID),
  cov   = as.numeric(cov),
  mthyl = as.numeric(mthyl)
)]
# if chrpos doesn't exist, create it
if (!("chrpos" %in% names(trail))) trail[, chrpos := paste0(chr, ":", pos)]
trail[, chrpos := as.character(chrpos)]

# 1) FAST subset to neuro CpGs
trail_sub <- trail[chrpos %chin% neuro_dt$chrpos]

# 2) Collapse duplicates -> 1 row per (chrpos, ID)
trail_per_id <- trail_sub[, .(
  cov  = sum(cov,  na.rm = TRUE),
  mthyl= sum(mthyl,na.rm = TRUE)
), by = .(chrpos, ID)]

trail_per_id[, beta := fifelse(cov > 0, mthyl / cov, NA_real_)]
trail_per_id <- trail_per_id[is.finite(beta)]

# ---- safe correlation helper ----
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

# 3) compute correlations for each cell
res_list <- vector("list", length(cov_grid) * length(app_grid))
k <- 1L

for (min_cov in cov_grid) {
  
  dt_cov <- trail_per_id[cov >= min_cov]
  
  # SD(beta) and number of IDs per CpG (after cov filter)
  stats_cov <- dt_cov[, .(
    n_ids   = .N,
    sd_beta = sd(beta, na.rm = TRUE)
  ), by = .(chrpos)]
  
  stats_cov <- stats_cov[is.finite(sd_beta)]
  
  for (min_app in app_grid) {
    
    keep <- stats_cov[n_ids >= min_app]
    ov <- merge(
      keep[, .(chrpos, sd_beta)],
      neuro_dt[, .(chrpos, SD_neuro)],
      by = "chrpos"
    )
    
    n_ov <- nrow(ov)
    ct <- safe_cor(ov$sd_beta, ov$SD_neuro, method = cor_method)
    
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

# order axes like your example
res[, min_cov := factor(min_cov, levels = cov_grid)]
res[, min_appear := factor(min_appear, levels = app_grid)]

# labels inside tiles (cor only, like the screenshot)
res[, lab := ifelse(is.na(cor), "NA", sprintf("%.2f", cor))]
# if you ALSO want counts in the tiles, swap to:
# res[, lab := ifelse(is.na(cor), "NA", sprintf("%.2f\n(n=%d)", cor, n_cpg))]

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
    x = "Minimum appearances (distinct IDs per CpG)",
    y = "Minimum coverage (per ID at CpG)"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())

print(p_heat)
ggsave("trail_neuro_sd_heatmap.png", p_heat, width = 8, height = 6, dpi = 300)

# optional: save the table behind the heatmap
fwrite(res, "trail_neuro_sd_heatmap_table.csv")












# ============================================================
# Heatmap of SPEARMAN correlation across:
#   min_cov (ID-level coverage)  x  min_appear (IDs per CpG)
# Comparing: Trail SD(beta) vs Neuro SD
# + shows n CpGs in each heatmap cell
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---- grid ----
cov_grid <- c(1, 2, 3, 5, 10)
app_grid <- c(2, 3, 5, 8, 10)
cor_method <- "spearman"

# ---- neuro prep ----
neuro_dt <- as.data.table(neuro)
neuro_dt[, CpG_Name := gsub(":", "_", as.character(CpG_Name))]
neuro_dt[, SD_neuro := as.numeric(Standard_Deviation)]
neuro_dt <- neuro_dt[is.finite(SD_neuro), .(CpG_Name, SD_neuro)]
neuro_dt[, chrpos := gsub("_", ":", CpG_Name)]
setkey(neuro_dt, chrpos)

# ---- trail prep ----
trail <- as.data.table(trail_petrous2)
trail[, `:=`(
  chr   = as.character(chr),
  pos   = as.integer(pos),
  ID    = as.character(ID),
  cov   = as.numeric(cov),
  mthyl = as.numeric(mthyl)
)]
if (!("chrpos" %in% names(trail))) trail[, chrpos := paste0(chr, ":", pos)]
trail[, chrpos := as.character(chrpos)]

# 1) FAST subset to neuro CpGs
trail_sub <- trail[chrpos %chin% neuro_dt$chrpos]

# 2) Collapse duplicates -> 1 row per (chrpos, ID)
trail_per_id <- trail_sub[, .(
  cov   = sum(cov, na.rm = TRUE),
  mthyl = sum(mthyl, na.rm = TRUE)
), by = .(chrpos, ID)]

trail_per_id[, beta := fifelse(cov > 0, mthyl / cov, NA_real_)]
trail_per_id <- trail_per_id[is.finite(beta)]

# ---- safe correlation helper ----
safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3) return(list(r = NA_real_, p = NA_real_))
  
  out <- tryCatch({
    ct <- cor.test(x, y, method = method)
    list(r = as.numeric(unname(ct$estimate)), p = as.numeric(ct$p.value))
  }, error = function(e) list(r = NA_real_, p = NA_real_))
  
  out
}

# 3) compute correlations for each cell
res_list <- vector("list", length(cov_grid) * length(app_grid))
k <- 1L

for (min_cov in cov_grid) {
  
  dt_cov <- trail_per_id[cov >= min_cov]
  
  # SD(beta) and number of IDs per CpG (after cov filter)
  stats_cov <- dt_cov[, .(
    n_ids   = .N,
    sd_beta = sd(beta, na.rm = TRUE)
  ), by = .(chrpos)]
  
  stats_cov <- stats_cov[is.finite(sd_beta)]
  
  for (min_app in app_grid) {
    
    keep <- stats_cov[n_ids >= min_app]
    
    ov <- merge(
      keep[, .(chrpos, sd_beta)],
      neuro_dt[, .(chrpos, SD_neuro)],
      by = "chrpos"
    )
    
    n_ov <- nrow(ov)
    ct <- safe_cor(ov$sd_beta, ov$SD_neuro, method = cor_method)
    
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

# label includes correlation + number of CpGs
res[, lab := ifelse(
  is.na(cor),
  sprintf("NA\n(n=%d)", n_cpg),
  sprintf("%.2f\n(n=%d)", cor, n_cpg)
)]

# ---- plot ----
p_heat <- ggplot(res, aes(x = min_appear, y = min_cov, fill = cor)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = lab), size = 3.6, lineheight = 0.95) +
  scale_fill_gradient2(
    name = "Spearman\nrho",
    midpoint = 0,
    low = "#0B1D39", mid = "#2C6AA3", high = "#5DB8FF",
    na.value = "grey90"
  ) +
  labs(
    title = "Method: Spearman",
    x = "number of inds (>=X)",
    y = "Minimum coverage (per ID at CpG)"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())

print(p_heat)
ggsave("trail_neuro_sd_heatmap_spearman.png", p_heat, width = 8, height = 6, dpi = 300)

# optional: save the table behind the heatmap
fwrite(res, "trail_neuro_sd_heatmap_table_spearman.csv")

# you should make a heatmap for spearman with the different coverages and stuff 
# do a plot which shows the variabiltiy between ancient and modern at the differnet coverages and stuff. 

