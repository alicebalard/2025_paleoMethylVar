## ============================================================
## Compare variability in Atlas10X (alpha) vs petrous methylation
## Using ALL CpGs (no Top/Bottom 100k filtering)
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
##   columns at least: name, alpha, chr, pos, pos2
## - cov_list_petrous: list of data.frames (one per petrous sample),
##   with columns:
##      V1 = chr (e.g. "chr1")
##      V2 = pos
##      V4 = methylation % (0–100)
##      V7 = coverage

Atlas_dt <- as.data.table(Atlas_dt)

## ------------------------------------------------------------------
## 1. Add chrpos to Atlas_dt (chr#:pos) for joining with petrous
## ------------------------------------------------------------------

# Re-derive chr and pos from 'name' to be safe, then chrpos
Atlas_dt[, c("chr_label", "pos_chr") := tstrsplit(name, "_", fixed = TRUE)]
Atlas_dt[, pos_chr := as.integer(pos_chr)]
Atlas_dt[, chrpos := paste0(chr_label, ":", pos_chr)]

## Quick sanity check
# head(Atlas_dt[, .(name, chr_label, pos_chr, chrpos, alpha)])

## ------------------------------------------------------------------
## 2. Build LONG petrous data for ALL CpGs
## ------------------------------------------------------------------

# You can adjust this coverage threshold if you like
min_cov_petrous <- 1

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
    cov >= min_cov_petrous
  )

petrous_long_dt <- as.data.table(petrous_long)

## ------------------------------------------------------------------
## 3. Per-CpG SD in each dataset (Atlas vs petrous)
## ------------------------------------------------------------------

## Atlas: SD of alpha across all Atlas observations per CpG
atlas_sd <- Atlas_dt[, .(
  n_obs_atlas = .N,
  sd_alpha    = sd(alpha, na.rm = TRUE),
  mean_alpha  = mean(alpha, na.rm = TRUE)
), by = chrpos]

## Petrous: SD of methylation across individuals per CpG
petrous_sd <- petrous_long_dt[, .(
  n_individuals = uniqueN(ID),
  sd_meth       = sd(meth, na.rm = TRUE),
  mean_meth     = mean(meth, na.rm = TRUE)
), by = chrpos]

## ------------------------------------------------------------------
## 4. Merge on shared CpGs and filter to CpGs with ≥2 obs in both
## ------------------------------------------------------------------

merged_sd <- merge(
  atlas_sd,
  petrous_sd,
  by  = "chrpos",
  all = FALSE
)

merged_sd_filt <- merged_sd %>%
  filter(
    n_obs_atlas   >= 2,
    n_individuals >= 2,
    !is.na(sd_alpha),
    !is.na(sd_meth)
  )

## Quick sanity check
# head(merged_sd_filt)
# summary(merged_sd_filt$sd_alpha)
# summary(merged_sd_filt$sd_meth)

## ------------------------------------------------------------------
## 5. Summaries: mean SD, SD of SD, SE of SD for each dataset
## ------------------------------------------------------------------

summary_long <- merged_sd_filt %>%
  transmute(
    chrpos,
    sd      = sd_alpha,
    dataset = "Atlas (alpha)"
  ) %>%
  bind_rows(
    merged_sd_filt %>%
      transmute(
        chrpos,
        sd      = sd_meth,
        dataset = "Petrous (meth %)"
      )
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

# y-limits based on 5th–95th percentile of all SDs
y_lims <- quantile(summary_long$sd, probs = c(0.05, 0.95), na.rm = TRUE)

## ------------------------------------------------------------------
## 6A. Plot: distribution of per-CpG SD (mean ± SE) for both datasets
## ------------------------------------------------------------------

p_dist <- ggplot() +
  # jittered per-CpG SDs
  geom_jitter(
    data = summary_long,
    aes(x = dataset, y = sd, colour = dataset),
    width = 0.15,
    alpha = 0.2,
    size  = 0.5,
    show.legend = FALSE
  ) +
  # mean ± SE per dataset
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
  geom_point(
    data = summary_sd,
    aes(x = dataset, y = mean_sd, colour = dataset),
    size = 3,
    show.legend = FALSE
  ) +
  coord_cartesian(ylim = y_lims) +
  labs(
    title = "Per-CpG variability: Atlas alpha vs petrous methylation",
    x = "",
    y = "SD (alpha or methylation %)"
  ) +
  theme_minimal(base_size = 13)

p_dist

## ------------------------------------------------------------------
## 6B. Scatter: SD_petrous vs SD_Atlas per CpG (shared CpGs)
## ------------------------------------------------------------------

p_scatter <- ggplot(merged_sd_filt,
                    aes(x = sd_meth, y = sd_alpha)) +
  geom_point(alpha = 0.3, size = 0.6) +
  labs(
    title = "Per-CpG SD: petrous methylation vs modern Atlas alpha",
    x = "SD of methylation (%) in petrous",
    y = "SD of alpha in Atlas"
  ) +
  theme_minimal(base_size = 13)

p_scatter
