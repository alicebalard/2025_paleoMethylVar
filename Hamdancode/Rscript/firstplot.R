suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(viridisLite)
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
    filter(!is.na(cov), cov >= min_cov) %>%
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

dist_all_cum <- map_df(1:10, ~cum_dist(trail_petrous, .x)) %>%
  mutate(label = factor(label, levels = paste0("≥", 1:10)))

# -------------------------------------------------------------------
# 2) Only methylated CpGs: "sites with methylation signal"
# -------------------------------------------------------------------
cum_dist_methylated <- function(dat, min_cov) {
  site_counts <- dat %>%
    filter(!is.na(cov), cov >= min_cov,
           !is.na(mthyl), mthyl > 0) %>%
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

# -------------------------------------------------------------------
# Dynamic y-axis range
# -------------------------------------------------------------------
y_max <- max(
  dist_all_cum$n_sites_at_least,
  dist_all_cum_meth$n_sites_at_least,
  na.rm = TRUE
)

y_break_max <- ceiling(log10(y_max))
y_breaks <- 10^(0:y_break_max)

# -------------------------------------------------------------------
# Clear paper-friendly palette
# -------------------------------------------------------------------
Hcolors <- setNames(
  viridis(10, option = "D", begin = 0.08, end = 0.95),
  paste0("≥", 1:10)
)

# -------------------------------------------------------------------
# Shared theme
# -------------------------------------------------------------------
paper_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    plot.title.position = "plot"
  )

# -------------------------------------------------------------------
# Plot 1: all CpGs
# -------------------------------------------------------------------
p1 <- ggplot(
  dist_all_cum,
  aes(x = n_people, y = n_sites_at_least, colour = label, group = label)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = Hcolors, name = "Min coverage") +
  scale_y_continuous(
    trans = "log10",
    breaks = y_breaks,
    labels = label_number(),
    limits = c(1, 10^y_break_max),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 10),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    title = "",
    x = "Sites shared by at least N individuals",
    y = "CpG sites (log10 scale)"
  ) +
  paper_theme

# -------------------------------------------------------------------
# Plot 2: methylated CpGs only
# -------------------------------------------------------------------
p2 <- ggplot(
  dist_all_cum_meth,
  aes(x = n_people, y = n_sites_at_least, colour = label, group = label)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = Hcolors, name = "Min coverage") +
  scale_y_continuous(
    trans = "log10",
    breaks = y_breaks,
    labels = label_number(),
    limits = c(1, 10^y_break_max),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 10),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    title = "",
    x = "Sites shared by at least N individuals",
    y = "Methylated CpG sites (log10 scale)"
  ) +
  paper_theme

# -------------------------------------------------------------------
# Combine side-by-side with one shared legend
# -------------------------------------------------------------------
combined_plot <- p1 + p2 +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# -------------------------------------------------------------------
# Save as high-resolution image for Word
# -------------------------------------------------------------------
ggsave(
  filename = "C:/Users/hamda/Desktop/UGI_thingy/petrous_cpg_combined.png",
  plot = combined_plot,
  width = 14,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "C:/Users/hamda/Desktop/UGI_thingy/petrous_cpg_combined.tiff",
  plot = combined_plot,
  width = 14,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)




library(dplyr)

site_summary <- trail_petrous_annotated %>%
  mutate(
    chr = as.character(chr),
    pos = as.integer(pos),
    mthyl = as.numeric(mthyl)
  ) %>%
  group_by(chr, pos) %>%
  summarise(
    any_methylation = any(!is.na(mthyl) & mthyl > 0),
    .groups = "drop"
  )

percent_sites_with_methylation <- mean(site_summary$any_methylation) * 100

cat("Percentage of sites with at least one methylation call:",
    round(percent_sites_with_methylation, 2), "%\n")