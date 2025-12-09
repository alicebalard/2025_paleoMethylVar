library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(scales)

#GET PROGRRESS BAR PACKAGE

#===============================================================================
# 0. Load CpG sets (top & bottom) and petrous coverage list
#===============================================================================

top100k    <- readRDS("3rdyear/genome/data/top100k_4Hamdan.rds")
bottom100k <- readRDS("3rdyear/genome/data/bottom100k_4Hamdan.rds")
# cov_list_petrous already in memory: list of data.frames

#===============================================================================
# 1. Add chrpos to top & bottom sets
#    (name looks like "chr1_819122" -> chr="chr1", pos=819122, chrpos="chr1:819122")
#===============================================================================

top100k <- top100k %>%
  mutate(
    chr    = sub("_(.*)", "", name),
    pos    = as.integer(sub(".*_", "", name)),
    chrpos = paste0(chr, ":", pos),
    group  = "Top 100k"
  )

bottom100k <- bottom100k %>%
  mutate(
    chr    = sub("_(.*)", "", name),
    pos    = as.integer(sub(".*_", "", name)),
    chrpos = paste0(chr, ":", pos),
    group  = "Bottom 100k"
  )

# All CpGs we care about (union of top & bottom sets)
cpg_union <- bind_rows(top100k, bottom100k) %>%
  distinct(chrpos, .keep_all = TRUE) %>%
  select(chrpos)

#===============================================================================
# 2. Long data frame from cov_list_petrous
#    - keep methylation %, coverage, and sample ID
#    - require coverage >= 3 (minimum depth)
#===============================================================================

petrous_long <- map_df(names(cov_list_petrous), function(id) {
  cov_list_petrous[[id]] %>%
    transmute(
      chr    = V1,
      pos    = V2,
      chrpos = paste0(V1, ":", V2),
      meth   = V4,   # percent methylation (0–100)
      cov    = V7,   # total coverage
      ID     = id    # sample/individual ID
    )
}) %>%
  # only CpGs that are in top or bottom 100k, and have depth >= 3
  filter(chrpos %in% cpg_union$chrpos,
         !is.na(cov), cov >= 2)

#===============================================================================
# 3. Find CpGs that are present (cov >= 3) in exactly 3 individuals
#===============================================================================

cpg_counts <- petrous_long %>%
  distinct(chrpos, ID) %>%          # one row per CpG × individual
  count(chrpos, name = "n_individuals")

cpg_3 <- cpg_counts %>%
  filter(n_individuals == 3) %>%    # present in exactly 3 samples
  select(chrpos)

#===============================================================================
# 4. Restrict to those CpGs and build LONG and WIDE tables
#===============================================================================

# LONG: one row per CpG × sample × (meth, cov)
cpg_3_long <- petrous_long %>%
  semi_join(cpg_3, by = "chrpos")

# Add top/bottom label + alpha (if you want it later)
cpg_3_long <- cpg_3_long %>%
  left_join(
    bind_rows(
      top100k    %>% select(chrpos, alpha, group),
      bottom100k %>% select(chrpos, alpha, group)
    ),
    by = "chrpos"
  )

# WIDE: rows = CpGs, columns = samples, values = methylation %
cpg_3_wide <- cpg_3_long %>%
  select(chrpos, ID, meth) %>%
  pivot_wider(names_from = ID, values_from = meth)

#===============================================================================
# 5. Per-CpG SD of methylation (and coverage) across those 3 individuals
#===============================================================================

cpg_3_sd <- cpg_3_long %>%
  group_by(chrpos, group) %>%   # group keeps Top vs Bottom separate
  summarise(
    n_individuals = n_distinct(ID),       # should be 3 for all
    sd_meth       = sd(meth, na.rm = TRUE),
    sd_cov        = sd(cov,  na.rm = TRUE),
    alpha         = first(alpha),         # alpha per CpG
    .groups       = "drop"
  )

# sanity check
head(cpg_3_long)
head(cpg_3_wide)
head(cpg_3_sd)
table(cpg_3_sd$n_individuals)  # should all be 3



library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(scales)

#===============================================================================
# 0. Prepare CpG sets: add chrpos and group label
#===============================================================================

top100k <- top100k %>%
  mutate(
    chr    = sub("_(.*)", "", name),
    pos    = as.integer(sub(".*_", "", name)),
    chrpos = paste0(chr, ":", pos),
    group  = "Top 100k"
  )

bottom100k <- bottom100k %>%
  mutate(
    chr    = sub("_(.*)", "", name),
    pos    = as.integer(sub(".*_", "", name)),
    chrpos = paste0(chr, ":", pos),
    group  = "Bottom 100k"
  )

#===============================================================================
# 1. Build long petrous coverage+methylation table
#    (only CpGs in top/bottom sets, and coverage >= 3)
#===============================================================================

# CpGs we care about (union of top & bottom)
cpg_keys <- bind_rows(top100k, bottom100k) %>%
  distinct(chrpos)

petrous_long <- map_df(names(cov_list_petrous), function(id) {
  cov_list_petrous[[id]] %>%
    transmute(
      chr    = V1,
      pos    = V2,
      chrpos = paste0(V1, ":", V2),
      meth   = V4,   # % methylation
      cov    = V7,   # coverage
      ID     = id    # individual/sample id
    )
}) %>%
  filter(
    chrpos %in% cpg_keys$chrpos,
    !is.na(cov),
    cov >= 1                # <<< coverage threshold (change here if needed)
  )

# attach Top/Bottom label
cpg_long_all <- petrous_long %>%
  left_join(
    bind_rows(
      top100k    %>% select(chrpos, group),
      bottom100k %>% select(chrpos, group)
    ),
    by = "chrpos"
  )

#===============================================================================
# 2. Per-CpG SD across individuals (for each group)
#===============================================================================

cpg_sd_all <- cpg_long_all %>%
  group_by(group, chrpos) %>%
  summarise(
    n_individuals = n_distinct(ID),
    sd_meth       = sd(meth, na.rm = TRUE),
    sd_cov        = sd(cov,  na.rm = TRUE),
    .groups       = "drop"
  )

# keep only CpGs present in 2 or more individuals
cpg_sd_plot <- cpg_sd_all %>%
  filter(n_individuals >= 2) %>%   # <-- changed here
  mutate(group = factor(group, levels = c("Bottom 100k", "Top 100k")))


#===============================================================================
# 3. Group summaries: mean SD, SD of SD, SE of SD
#===============================================================================

summary_sd <- cpg_sd_plot %>%
  group_by(group) %>%
  summarise(
    mean_sd = mean(sd_meth, na.rm = TRUE),
    sd_sd   = sd(sd_meth,   na.rm = TRUE),      # total SD of CpG SDs
    n_cpgs  = n(),
    se_sd   = sd_sd / sqrt(n_cpgs),             # SE of mean SD
    .groups = "drop"
  )

print(summary_sd)

# nice y-limits for zooming if you want them (5th–95th percentile)
y_lims <- quantile(cpg_sd_plot$sd_meth, probs = c(0.05, 0.95), na.rm = TRUE)

#===============================================================================
# 4A. Plot 1: mean SD ± SE (uncertainty of the mean)
#===============================================================================

p_se <- ggplot() +
  # jittered individual CpG SDs
  geom_jitter(
    data = cpg_sd_plot,
    aes(x = group, y = sd_meth, colour = group),
    width = 0.15,
    alpha = 0.25,
    size  = 0.8,
    show.legend = FALSE
  ) +
  # mean ± SE
  geom_errorbar(
    data = summary_sd,
    aes(x = group,
        ymin = mean_sd - se_sd,
        ymax = mean_sd + se_sd,
        colour = group),
    width = 0.12,
    size  = 0.9,
    show.legend = FALSE
  ) +
  geom_point(
    data = summary_sd,
    aes(x = group, y = mean_sd, colour = group),
    size = 3,
    show.legend = FALSE
  ) +
  coord_cartesian(ylim = y_lims) +
  labs(
    title = "SD of methylation across individuals\nTop vs Bottom 100k CpGs (mean ± SE)",
    x = "CpG set",
    y = "SD of methylation (%)"
  ) +
  theme_minimal(base_size = 13)

p_se

#===============================================================================
# 4B. Plot 2: mean SD ± SD (spread of CpG SDs)
#===============================================================================

p_sd <- ggplot() +
  geom_jitter(
    data = cpg_sd_plot,
    aes(x = group, y = sd_meth, colour = group),
    width = 0.15,
    alpha = 0.25,
    size  = 0.8,
    show.legend = FALSE
  ) +
  # mean ± total SD
  geom_errorbar(
    data = summary_sd,
    aes(x = group,
        ymin = pmax(mean_sd - sd_sd, 0),  # avoid going below 0
        ymax = mean_sd + sd_sd,
        colour = group),
    width = 0.12,
    size  = 0.9,
    show.legend = FALSE
  ) +
  geom_point(
    data = summary_sd,
    aes(x = group, y = mean_sd, colour = group),
    size = 3,
    show.legend = FALSE
  ) +
  labs(
    title = "SD of methylation across individuals\nTop vs Bottom 100k CpGs (mean ± SD)",
    x = "CpG set",
    y = "SD of methylation (%)"
  ) +
  theme_minimal(base_size = 13)

p_sd




#===============================================================================
# 4B. Plot 2: mean SD ± SD (spread of CpG SDs)
#===============================================================================

p_sd2 <- ggplot() +
  # mean ± total SD
  geom_errorbar(
    data = summary_sd,
    aes(x = group,
        ymin = pmax(mean_sd - sd_sd, 0),  # avoid going below 0
        ymax = mean_sd + sd_sd,
        colour = group),
    width = 0.12,
    size  = 0.9,
    show.legend = FALSE
  ) +
  geom_point(
    data = summary_sd,
    aes(x = group, y = mean_sd, colour = group),
    size = 3,
    show.legend = FALSE
  ) +
  labs(
    title = "SD of methylation across individuals\nTop vs Bottom 100k CpGs (mean ± SD)",
    x = "CpG set",
    y = "SD of methylation (%)"
  ) +
  theme_minimal(base_size = 13)

p_sd2
