library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)

## 1. hg38 CpGs we care about -------------------------------------

cpgs_hg38 <- Mdata_hg38 %>%
  filter(!is.na(chr_hg38), !is.na(hg38.pos)) %>%
  mutate(chrpos = paste0(chr_hg38, ":", hg38.pos)) %>%
  distinct(chrpos)

## 2. Long DF of coverage, restricted to Mdata_hg38 CpGs ----------

petrous_long <- map_df(names(cov_list_petrous), function(id) {
  cov_list_petrous[[id]] %>%
    transmute(
      chr    = V1,
      pos    = V2,
      chrpos = paste0(V1, ":", V2),
      cov    = V7,
      ID     = id
    )
})

petrous_hg38 <- petrous_long %>%
  semi_join(cpgs_hg38, by = "chrpos")

## 3. Cumulative "≥ N people" function (safe) ---------------------

cum_dist <- function(dat, min_cov) {
  site_counts <- dat %>%
    filter(!is.na(cov), cov >= min_cov) %>%
    distinct(chrpos, ID) %>%
    count(chrpos, name = "n_people")
  
  if (nrow(site_counts) == 0L) {
    return(tibble(
      n_people    = integer(),
      n_sites     = integer(),
      n_sites_cum = integer()
    ))
  }
  
  site_counts %>%
    count(n_people, name = "n_sites") %>%
    arrange(n_people) %>%
    complete(
      n_people = full_seq(n_people, 1),
      fill = list(n_sites = 0)
    ) %>%
    mutate(
      n_sites_cum = rev(cumsum(rev(n_sites)))
    )
}

## 4. Choose min-coverage thresholds that make sense --------------

# Look at max coverage to decide how far to go
max_cov <- max(petrous_hg38$cov, na.rm = TRUE)
max_cov

# Here I’ll use 1:3; you can bump up to max_cov if there’s signal
min_cov_vec <- 1:7

dist_long <- map_df(min_cov_vec, function(mc) {
  cum_dist(petrous_hg38, mc) %>%
    mutate(min_cov = paste0("≥", mc))
})

dist_long$min_cov <- factor(
  dist_long$min_cov,
  levels = paste0("≥", min_cov_vec)
)

## 5. Plot in your original style ---------------------------------

ggplot(dist_long,
       aes(x = n_people, y = n_sites_cum, fill = min_cov)) +
  geom_col(colour = "black", width = 1, position = "identity") +
  scale_y_log10(labels = label_comma()) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Mdata_hg38 CpG sites shared by at least N people\nat different coverage thresholds",
    x = "N people (threshold)",
    y = "CpG sites (log scale)",
    fill = "Min coverage"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right"
  )
