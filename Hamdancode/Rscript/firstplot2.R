suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(scales)
})

# ---------------------------
# Prep: ensure correct types
# ---------------------------
trail_petrous <- trail_petrous_annotated %>%
  mutate(
    chr   = as.character(chr),
    pos   = as.integer(pos),
    ID    = as.character(ID),
    cov   = as.numeric(cov),
    mthyl = as.numeric(mthyl)
  )

# -------------------------------------------------------------------
# 1) All CpGs: "sites present at coverage ≥ X in at least N people"
# -------------------------------------------------------------------
cum_dist <- function(dat, min_cov) {
  site_counts <- dat %>%
    filter(!is.na(cov), cov >= min_cov) %>%    # cumulative: ≥ min_cov
    distinct(chr, pos, ID) %>%                 # one row per site/person
    count(chr, pos, name = "n_people")         # people per site
  
  dist <- site_counts %>%
    count(n_people, name = "n_sites") %>%      # sites with exactly N people
    arrange(n_people) %>%
    complete(n_people = full_seq(n_people, 1), fill = list(n_sites = 0)) %>%
    mutate(
      n_sites_at_least = rev(cumsum(rev(n_sites))),  # sites with ≥ N people
      n_people         = as.integer(n_people),
      min_cov          = min_cov,
      label            = paste0("≥", min_cov)
    )
  
  dist %>% select(n_people, n_sites_at_least, min_cov, label)
}

dist_all_cum <- map_df(1:10, ~cum_dist(trail_petrous, .x)) %>%
  mutate(label = factor(label, levels = paste0("≥", 1:10)))

Hcolors <- c(
  "≥1"  = "#F2F2F2",
  "≥2"  = "#E0E0E0",
  "≥3"  = "#CCCCCC",
  "≥4"  = "#B3B3B3",
  "≥5"  = "#999999",
  "≥6"  = "#808080",
  "≥7"  = "#666666",
  "≥8"  = "#4D4D4D",
  "≥9"  = "#333333",
  "≥10" = "#1A1A1A"
)

pdf("C:/Users/hamda/Desktop/UGI_thingy/fig1_petrous_cov1.pdf", width = 10, height = 5)
ggplot() +
  geom_col(
    data = dist_all_cum,
    aes(x = n_people, y = n_sites_at_least, fill = label),
    position = "identity", width = 0.9, alpha = 0.65, colour = "black"
  ) +
  scale_fill_manual(values = Hcolors, name = "Min coverage") +
  scale_y_continuous(
    trans = "log10",
    labels = label_number(scale_cut = cut_short_scale()),
    breaks = 10^(0:6),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 12),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(
    title = "Petrous CpG sites shared by at least N people at different coverage thresholds",
    x = "N people (threshold)",
    y = "CpG sites (log scale)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    legend.position = "right"
  )
dev.off()

# -------------------------------------------------------------------
# 2) Only methylated CpGs: "sites with methylation signal"
# -------------------------------------------------------------------
cum_dist_methylated <- function(dat, min_cov) {
  site_counts <- dat %>%
    filter(!is.na(cov), cov >= min_cov,
           !is.na(mthyl), mthyl > 0) %>%       # methylated only
    distinct(chr, pos, ID) %>%
    count(chr, pos, name = "n_people")
  
  dist <- site_counts %>%
    count(n_people, name = "n_sites") %>%
    arrange(n_people) %>%
    complete(n_people = full_seq(n_people, 1), fill = list(n_sites = 0)) %>%
    mutate(
      n_sites_at_least = rev(cumsum(rev(n_sites))),
      n_people         = as.integer(n_people),
      min_cov          = min_cov,
      label            = paste0("≥", min_cov)
    )
  
  dist %>% select(n_people, n_sites_at_least, min_cov, label)
}

dist_all_cum_meth <- map_df(1:10, ~cum_dist_methylated(trail_petrous, .x)) %>%
  mutate(label = factor(label, levels = paste0("≥", 1:10)))

pdf("C:/Users/hamda/Desktop/UGI_thingy/fig2_petrous_methylated_cov1.pdf", width = 10, height = 5)
ggplot() +
  geom_col(
    data = dist_all_cum_meth,
    aes(x = n_people, y = n_sites_at_least, fill = label),
    position = "identity", width = 0.9, alpha = 0.65, colour = "black"
  ) +
  scale_fill_manual(values = Hcolors, name = "Min coverage") +
  scale_y_continuous(
    trans = "log10",
    labels = label_number(scale_cut = cut_short_scale()),
    breaks = 10^(0:6),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 12),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(
    title = "Petrous methylated CpG sites shared by at least N people (by coverage threshold)",
    x = "N people",
    y = "Methylated CpG sites (log scale)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    legend.position = "right"
  )
dev.off()