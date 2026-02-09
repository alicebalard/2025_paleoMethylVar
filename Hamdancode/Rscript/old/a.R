library(dplyr)
library(tidyr)
library(ggplot2)
head(testmerge)

N3testmerge <- table(testmerge$V1)
N3testmerge <- as.data.frame(N3testmerge)
N3testmerge <- subset(N3testmerge, Freq >= 5)

testmerge101 <- testmerge[testmerge$V1 %in% N3testmerge$Var1, ]
head(testmerge101)
#-- 1a. Split V1 so Chr and Pos are separate
nrow(testmerge101)
nrow(testmerge)

testmerge2 <- testmerge101 %>%
  separate(V1, into = c("Chr", "Pos"), sep = "_", convert = TRUE) %>%
  mutate(SiteID = paste0("chr", Chr, "_", Pos))   # convenient single key




site_var <- testmerge2 %>%
  group_by(SiteID, period4ana) %>%
  filter(n() > 1) %>%                           # keep sites seen ≥ 2× in that period
  summarise(var_site = var(V5, na.rm = TRUE),
            n_ind    = n(),
            .groups  = "drop")



ggplot(site_var, aes(x = var_site)) +
  geom_histogram(binwidth = 0.01, boundary = 0, fill = "#69b3a2") +
  facet_wrap(~ period4ana, scales = "free_y", nrow = 1) +
  scale_x_continuous(
    breaks = seq(0, 0.30, by = 0.1),   # ticks at 0.00, 0.01, 0.02, …
    limits = c(0, 0.30),                # zoom if you only care up to 0.30
    expand = expansion(mult = 0)        # trim the empty margins
  ) +
  labs(x = "Site-specific variance of methylation (V5)",
       y = "Number of CpG sites") +
  theme_minimal()
      


site_var$var_site[site_var$var_site == 0] <- 1e-5
ggplot(site_var, aes(x = period4ana, y = var_site)) +
  geom_boxplot(outlier.size = 0.5, fill = "#2ca25f") +
  scale_y_continuous(trans = "log10") +          # optional log-scale
  labs(x = "Period", y = "Site-specific variance (log10)") +
  theme_minimal()


site_var_wide <- site_var %>%
  pivot_wider(names_from  = period4ana,
              values_from = var_site)

site_var_wide


site_changes <- site_var_wide %>%
  mutate(delta = Imperial - `Iron/Republic`) %>%  # change sign/period names as needed
  filter(!is.na(delta)) %>%
  arrange(desc(delta))





head(site_changes, 20)      # top 20 sites: variance went up the most
tail(site_changes, 20) 




top_delta <- site_changes %>% slice_head(n = 8000)   # pick a threshold you like
top_delta
ggplot(top_delta,
       aes(x = reorder(SiteID, delta), y = delta)) +
  geom_col(fill = "#d95f02") +
  coord_flip() +
  labs(x = "CpG site",
       y = "Δ variance (Imperial minus Iron/Republic)") +
  theme_minimal()



library(purrr)

# choose two periods to compare
period_pair <- c("Iron/Republic", "Imperial")

var_test_df <- testmerge2 %>%
  filter(period4ana %in% period_pair) %>%
  group_by(SiteID) %>%
  nest() %>%
  mutate(p_val = map_dbl(data, ~ {
    d <- .
    v1 <- d$V5[d$period4ana == period_pair[1]]
    v2 <- d$V5[d$period4ana == period_pair[2]]
    if(length(v1) > 1 && length(v2) > 1) {
      var.test(v1, v2)$p.value     # F-test
    } else NA_real_
  })) %>%
  ungroup() %>%
  filter(!is.na(p_val)) %>%
  mutate(p_adj = p.adjust(p_val, method = "BH"))    # multiple-testing correction



t1t <- filter(var_test_df, p_adj < 0.05) %>%
  arrange(p_adj) %>%
  select(SiteID, p_adj) %>%
  head()
 
t1t

t1t_full <- var_test_df %>%
  filter(SiteID %in% t1t$SiteID) %>%
  mutate(p_val_fmt  = formatC(p_val,  format = "e", digits = 2),
         p_adj_fmt  = formatC(p_adj,  format = "e", digits = 2))

t1t_full %>% 
  select(SiteID, p_val_fmt, p_adj_fmt)

# pull only the six significant sites
sig_data <- testmerge2 %>%
  filter(SiteID %in% t1t$SiteID,
         period4ana %in% period_pair)

# summarise variance (and mean) in each period
sig_summary <- sig_data %>%
  group_by(SiteID, period4ana) %>%
  summarise(var  = var(V5, na.rm = TRUE),
            mean = mean(V5, na.rm = TRUE),
            n    = n(),
            .groups = "drop") %>%
  pivot_wider(names_from = period4ana,
              values_from = c(var, mean, n),
              names_sep   = ".")

sig_summary
ggplot(sig_data, aes(x = period4ana, y = V5, colour = period4ana)) +
  geom_jitter(width = 0.15, height = 0, size = 1.2, alpha = 0.7) +
  facet_wrap(~ SiteID, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = "Methylation β (V5)",
       title = "CpG sites with significantly different variance") +
  theme_minimal() +
  theme(legend.position = "none")






library(tidyverse)

# 1. Pivot to wide format
site_var_wide <- site_var %>%
  pivot_wider(names_from  = period4ana,
              values_from = var_site)

# 2. Get all pairwise combinations of periods
periods <- setdiff(names(site_var_wide), "SiteID")
period_pairs <- combn(periods, 2, simplify = FALSE)

# 3. Loop through each pair and plot top 50 differences
for (pair in period_pairs) {
  period1 <- pair[1]
  period2 <- pair[2]
  
  comp_name <- paste0(period1, "_vs_", period2)
  
  site_changes <- site_var_wide %>%
    mutate(delta = .data[[period1]] - .data[[period2]]) %>%
    filter(!is.na(delta)) %>%
    arrange(desc(delta))
  
  top_delta <- site_changes %>% slice_max(delta, n = 50)
  
  p <- ggplot(top_delta,
              aes(x = reorder(SiteID, delta), y = delta)) +
    geom_col(fill = "#d95f02") +
    coord_flip() +
    labs(
      title = paste("Top Δ variance:", period1, "minus", period2),
      x = "CpG site",
      y = paste("Δ variance (", period1, "minus", period2, ")")
    ) +
    theme_minimal()
  
  print(p)
}

# 1. Pivot to wide format
site_var_wide <- site_var %>%
  pivot_wider(names_from  = period4ana,
              values_from = var_site)

# 2. Get all pairwise combinations of periods
periods <- setdiff(names(site_var_wide), "SiteID")
period_pairs <- combn(periods, 2, simplify = FALSE)

# 3. Loop through each pair and plot top 50 differences
for (pair in period_pairs) {
  period1 <- pair[1]
  period2 <- pair[2]
  
  comp_name <- paste0(period1, "_vs_", period2)
  
  site_changes <- site_var_wide %>%
    filter(.data[[period1]] != 0, .data[[period2]] != 0) %>%
    mutate(delta = .data[[period1]] - .data[[period2]]) %>%
    filter(!is.na(delta)) %>%
    arrange(desc(delta))
  
  top_delta <- site_changes %>% slice_max(delta, n = 1000)
  
  p <- ggplot(top_delta,
              aes(x = reorder(SiteID, delta), y = delta)) +
    geom_col(fill = "#d95f02") +
    coord_flip() +
    labs(
      title = paste("Top Δ variance:", period1, "minus", period2),
      x = "CpG site",
      y = paste("Δ variance (", period1, "minus", period2, ")")
    ) +
    theme_minimal()
  
  print(p)
}




