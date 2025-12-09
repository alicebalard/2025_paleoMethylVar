install.packages("viridis")  # run once
library(viridis)


library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(scales)
library(viridis)

out_dir <- "C:/Users/hamda/Desktop/UGI_thingy"

#-----------------------------------------------------------
# 3A. Core cumulative distribution function
#-----------------------------------------------------------

cum_dist <- function(dat, min_cov) {
  site_counts <- dat %>%
    filter(!is.na(cov), cov >= min_cov) %>%           # ≥ min_cov
    distinct(chr, pos, ID) %>%                        # one row per site/person
    count(chr, pos, name = "n_people")                # people per CpG
  
  dist <- site_counts %>%
    count(n_people, name = "n_sites") %>%             # CpGs with exactly N people
    arrange(n_people) %>%
    complete(n_people = full_seq(n_people, 1),
             fill = list(n_sites = 0)) %>%
    mutate(
      n_sites_at_least = rev(cumsum(rev(n_sites))),   # CpGs with ≥ N people
      n_people         = as.integer(n_people),
      min_cov          = min_cov,
      label            = paste0("≥", min_cov)
    )
  
  dist %>% select(n_people, n_sites_at_least, min_cov, label)
}

make_cum_data <- function(trail_df, max_cov = 10) {
  map_df(1:max_cov, ~cum_dist(trail_df, .x)) %>%
    mutate(label = factor(label, levels = paste0("≥", 1:max_cov)))
}

#-----------------------------------------------------------
# 3B. Generic plotting function
#-----------------------------------------------------------

make_cum_plot <- function(cum_df, title, filename) {
  p <- ggplot(cum_df,
              aes(x = n_people,
                  y = n_sites_at_least,
                  fill = label)) +
    geom_col(position = "identity", width = 0.9, alpha = 0.9) +
    scale_fill_viridis_d(
      option   = "C",
      direction = -1,
      name     = "Min coverage"
    ) +
    scale_y_continuous(
      trans   = "log10",
      labels  = label_number(scale_cut = cut_short_scale()),
      breaks  = 10^(0:6),
      expand  = expansion(mult = c(0.02, 0.08))
    ) +
    scale_x_continuous(
      breaks = pretty_breaks(n = 12),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title = title,
      x     = "N people (threshold)",
      y     = "CpG sites (log scale)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor   = element_blank(),
      plot.title.position = "plot",
      legend.position     = "right"
    )
  
  print(p)
  
  ggsave(
    filename = file.path(out_dir, filename),
    plot     = p,
    width    = 10,
    height   = 6
  )
}

#-----------------------------------------------------------
# 3C. Compute and plot for:
#     - molar-only dataset
#     - petrous-only dataset
#-----------------------------------------------------------

dist_molar   <- make_cum_data(trail_molar,   max_cov = 10)
dist_petrous <- make_cum_data(trail_petrous, max_cov = 10)

make_cum_plot(
  dist_molar,
  title    = "CpG sites shared by at least N people (Vac molar only + all other samples)",
  filename = "Fig_cum_molar_only.pdf"
)

make_cum_plot(
  dist_petrous,
  title    = "CpG sites shared by at least N people (Vac petrous only + all other samples)",
  filename = "Fig_cum_petrous_only.pdf"
)

#-----------------------------------------------------------
# (Optional) versions restricted to CpGs with methylation signal
#-----------------------------------------------------------

cum_dist_methylated <- function(dat, min_cov) {
  site_counts <- dat %>%
    filter(!is.na(cov), cov >= min_cov,
           !is.na(mthyl), mthyl > 0) %>%        # require methylation signal
    distinct(chr, pos, ID) %>%
    count(chr, pos, name = "n_people")
  
  dist <- site_counts %>%
    count(n_people, name = "n_sites") %>%
    arrange(n_people) %>%
    complete(n_people = full_seq(n_people, 1),
             fill = list(n_sites = 0)) %>%
    mutate(
      n_sites_at_least = rev(cumsum(rev(n_sites))),
      n_people         = as.integer(n_people),
      min_cov          = min_cov,
      label            = paste0("≥", min_cov)
    )
  
  dist %>% select(n_people, n_sites_at_least, min_cov, label)
}

make_cum_data_meth <- function(trail_df, max_cov = 10) {
  map_df(1:max_cov, ~cum_dist_methylated(trail_df, .x)) %>%
    mutate(label = factor(label, levels = paste0("≥", 1:max_cov)))
}

dist_molar_meth   <- make_cum_data_meth(trail_molar,   max_cov = 10)
dist_petrous_meth <- make_cum_data_meth(trail_petrous, max_cov = 10)

make_cum_plot(
  dist_molar_meth,
  title    = "Methylated CpG sites shared by at least N people (Vac molar only + others)",
  filename = "Fig_cum_molar_only_methylated.pdf"
)

make_cum_plot(
  dist_petrous_meth,
  title    = "Methylated CpG sites shared by at least N people (Vac petrous only + others)",
  filename = "Fig_cum_petrous_only_methylated.pdf"
)
