## ============================================================
## Compare variability in Atlas10X (alpha) vs petrous methylation
## Using ONLY CpGs present in BOTH datasets
## ============================================================

library(data.table)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(scales)

## ------------------------------------------------------------------
## 0. Inputs (assumed already loaded)
## ------------------------------------------------------------------
## - Atlas_dt from prepAtlasdt()
##   columns at least: name, alpha, chr, pos, pos2  (one row per CpG)
## - cov_list_petrous: list of data.frames (one per petrous sample),
##   with columns:
##      V1 = chr (e.g. "chr1")
##      V2 = pos
##      V4 = methylation % (0–100)
##      V7 = coverage

Atlas_dt <- as.data.table(Atlas_dt)

## ------------------------------------------------------------------
## 1. Make chrpos in BOTH datasets and find shared CpGs
## ------------------------------------------------------------------

## 1A. chrpos in Atlas
## If prepAtlasdt() already gave you chr and pos, you could just:
## Atlas_dt[, chrpos := paste0(chr, ":", pos)]
## Below is a safe version using 'name', which should look like "chr1_819122"

Atlas_dt[, c("chr_label", "pos_chr") := tstrsplit(name, "_", fixed = TRUE)]
Atlas_dt[, pos_chr := as.integer(pos_chr)]
Atlas_dt[, chrpos := paste0(chr_label, ":", pos_chr)]

## 1B. Build LONG petrous data (all CpGs) with chrpos

petrous_long <- map_df(names(cov_list_petrous), function(id) {
  cov_list_petrous[[id]] %>%
    transmute(
      chr    = V1,
      pos    = V2,
      chrpos = paste0(V1, ":", V2),
      meth   = V4,   # % methylation
      cov    = V7,   # coverage
      ID     = id    # individual/sample ID
    )
}) %>%
  filter(
    !is.na(cov),
    cov >= 3
  )

petrous_long_dt <- as.data.table(petrous_long)

## 1C. Find CpGs present in BOTH Atlas and petrous
atlas_chrpos   <- unique(Atlas_dt$chrpos)
petrous_chrpos <- unique(petrous_long_dt$chrpos)

shared_chrpos <- intersect(atlas_chrpos, petrous_chrpos)
length(shared_chrpos)  # how many CpGs in common?

## ------------------------------------------------------------------
## 2. SUBSET both datasets to SHARED CpGs only
## ------------------------------------------------------------------

Atlas_shared   <- Atlas_dt[chrpos %in% shared_chrpos]
petrous_shared <- petrous_long_dt[chrpos %in% shared_chrpos]

## Quick sanity check
length(unique(Atlas_shared$chrpos))
length(unique(petrous_shared$chrpos))





## ------------------------------------------------------------------
## 2.5 NEW FILTER (requested):
## Keep only CpGs that have >=4 occurrences in petrous_long (after cov>=3)
## Here, "occurrences" = number of rows/measurements for that CpG in petrous_long,
## which in your setup is effectively "number of individuals with data".
## ------------------------------------------------------------------

petrous_counts <- petrous_shared[, .(n_occ = .N), by = chrpos]
keep_chrpos_4  <- petrous_counts[n_occ >= 4, chrpos]

petrous_shared <- petrous_shared[chrpos %in% keep_chrpos_4]
Atlas_shared   <- Atlas_shared[chrpos %in% keep_chrpos_4]

## sanity checks after the >=4 filter
cat("CpGs after shared filter:", length(shared_chrpos), "\n")
cat("CpGs after petrous>=4 occurrences filter:", length(unique(petrous_shared$chrpos)), "\n")



## ------------------------------------------------------------------
## 3. Per-CpG SD in each dataset (on SHARED CpGs only)
## ------------------------------------------------------------------

## 3A. Atlas: SD of alpha across observations per CpG
## NOTE: with your current Atlas_dt, there is only ONE row per CpG,
## so n_obs_atlas will be 1 and sd_alpha will be NA.
## We still compute it in case you later add replicates;
## mean_alpha is the main hv summary right now.

atlas_sd <- Atlas_shared[, .(
  n_obs_atlas = .N,
  sd_alpha    = sd(alpha, na.rm = TRUE),
  mean_alpha  = mean(alpha, na.rm = TRUE)
), by = chrpos]

## 3B. Petrous: SD of methylation across individuals per CpG
petrous_sd <- petrous_shared[, .(
  n_individuals = uniqueN(ID),
  sd_meth       = sd(meth, na.rm = TRUE),
  mean_meth     = mean(meth, na.rm = TRUE)
), by = chrpos]

## ------------------------------------------------------------------
## 4. Merge per-CpG summaries and FILTER
## ------------------------------------------------------------------

merged_sd <- merge(
  atlas_sd,
  petrous_sd,
  by  = "chrpos",
  all = FALSE
)

## Here we ONLY require >=2 petrous individuals, because Atlas has 1 obs/CpG.
## If you ever have multiple Atlas observations per CpG, you can add
## n_obs_atlas >= 2 and !is.na(sd_alpha) back in.
merged_sd_filt <- merged_sd %>%
  filter(
    n_individuals >= 2,
    !is.na(sd_meth)
    # sd_alpha will be NA with 1 obs; keep mean_alpha instead
  )

## Quick sanity check
# head(merged_sd_filt)
# summary(merged_sd_filt$mean_alpha)
# summary(merged_sd_filt$sd_meth)

## ------------------------------------------------------------------
## 5. Summaries: mean SD, SD of SD, SE of SD for each dataset
##    (using SD of alpha only if it exists; otherwise petrous side)
## ------------------------------------------------------------------

## For now, the only meaningful SDs are on the petrous side;
## but we keep the same structure you had.

summary_long <- merged_sd_filt %>%
  transmute(
    chrpos,
    sd      = sd_meth,
    dataset = "Petrous (meth %)"
  )

summary_sd <- summary_long %>%
  group_by(dataset) %>%
  summarise(
    mean_sd = mean(sd, na.rm = TRUE),
    sd_sd   = sd(sd,   na.rm = TRUE),
    n_cpgs  = n(),
    se_sd   = sd_sd / sqrt(n_cpgs),
    .groups = "drop"
  )

print(summary_sd)

# y-limits based on 5th–95th percentile of SDs (petrous only for now)
y_lims <- quantile(summary_long$sd, probs = c(0.05, 0.95), na.rm = TRUE)

## ------------------------------------------------------------------
## 6A. Plot: distribution of per-CpG SD(meth) (mean ± SE)
##         (error bars only)
## ------------------------------------------------------------------

p_dist <- ggplot() +
  geom_errorbar(
    data = summary_sd,
    aes(x = dataset,
        ymin = mean_sd - se_sd,
        ymax = mean_sd + se_sd,
        colour = dataset),
    width = 0.12,
    size  = 0.9,
    show.legend = FALSE
  ) +
  coord_cartesian(ylim = y_lims) +
  labs(
    title = "Per-CpG variability in petrous methylation\n(shared CpGs with Atlas)",
    x = "",
    y = "SD of methylation (%)"
  ) +
  theme_minimal(base_size = 13)

p_dist






p_dist <- ggplot() +
  geom_errorbar(
    data = summary_sd,
    aes(
      x    = dataset,
      ymin = mean_sd - sd_sd,   # <- use SD here
      ymax = mean_sd + sd_sd,   # <- and here
      colour = dataset
    ),
    width = 0.12,
    size  = 0.9,
    show.legend = FALSE
  ) +
  coord_cartesian(ylim = y_lims) +
  labs(
    title = "Per-CpG variability in petrous methylation\n(shared CpGs with Atlas)",
    x = "",
    y = "SD of methylation (%)"   # optionally clarify: "mean ± 1 SD"
  ) +
  theme_minimal(base_size = 13)




p_dist

## ------------------------------------------------------------------
## 6B. Scatter: SD_petrous vs Atlas alpha for shared CpGs
##         (regression line only)
## ------------------------------------------------------------------

p_scatter <- ggplot(merged_sd_filt,
                    aes(x = sd_meth, y = mean_alpha)) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Petrous inter-individual variability vs modern Atlas alpha\n(shared CpGs only)",
    x = "SD of methylation (%) in petrous",
    y = "Atlas alpha (mean per CpG)"
  ) +
  theme_minimal(base_size = 13)

p_scatter




p_scatter2 <- ggplot(merged_sd_filt,
                    aes(x = sd_meth, y = mean_alpha)) +
  geom_point(alpha = 0.3, size = 0.7) +      # one point per CpG
  geom_smooth(method = "lm", se = FALSE) +   # regression line only
  labs(
    title = "Petrous inter-individual variability vs modern Atlas alpha\n(shared CpGs only)",
    x = "SD of methylation (%) in petrous",
    y = "Atlas alpha (mean per CpG)"
  ) +
  theme_minimal(base_size = 13)

p_scatter2








