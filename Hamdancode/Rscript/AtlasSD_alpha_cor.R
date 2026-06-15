suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

dir.create("plots", showWarnings = FALSE)

# --- checks ---
stopifnot(exists("AtlasSD"), exists("Atlas_annotated"))

# --- prep AtlasSD ---
AtlasSD_dt <- as.data.table(AtlasSD)
stopifnot(all(c("CpG_Name","Standard_Deviation") %in% names(AtlasSD_dt)))

AtlasSD_dt[, CpG_Name := as.character(CpG_Name)]
AtlasSD_dt[, SD_Atlas := as.numeric(Standard_Deviation)]

# if duplicates exist, average SD per CpG
AtlasSD_dt <- AtlasSD_dt[is.finite(SD_Atlas),
                         .(SD_Atlas = mean(SD_Atlas, na.rm = TRUE)),
                         by = CpG_Name]
setkey(AtlasSD_dt, CpG_Name)

# --- prep Atlas_annotated (alpha) ---
Ann_dt <- as.data.table(Atlas_annotated)
stopifnot(all(c("name","alpha") %in% names(Ann_dt)))

Ann_dt[, CpG_Name := as.character(name)]
Ann_dt[, alpha := as.numeric(alpha)]

# if duplicates exist, average alpha per CpG
Ann_dt <- Ann_dt[is.finite(alpha) & !is.na(CpG_Name),
                 .(alpha = mean(alpha, na.rm = TRUE)),
                 by = CpG_Name]
setkey(Ann_dt, CpG_Name)

# --- overlap and correlation ---
ov <- merge(Ann_dt, AtlasSD_dt, by = "CpG_Name")
cat("Overlapping CpGs:", nrow(ov), "\n")

r <- NA_real_
p <- NA_real_
if (nrow(ov) >= 3) {
  ct <- suppressWarnings(cor.test(ov$alpha, ov$SD_Atlas, method = "pearson"))
  r <- unname(ct$estimate)
  p <- ct$p.value
}

# --- LINE ONLY plot (no points) ---
p_line <- ggplot(ov, aes(x = alpha, y = SD_Atlas)) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  labs(
    title = "Atlas SD vs Atlas alpha (overlapping CpGs)",
    subtitle = if (nrow(ov) >= 3)
      sprintf("n=%d | Pearson r=%.3f | p=%.3g", nrow(ov), r, p)
    else
      sprintf("n=%d | Pearson r=NA (need n>=3)", nrow(ov)),
    x = "alpha (Atlas_annotated)",
    y = "Standard_Deviation (AtlasSD)"
  ) +
  theme_minimal(base_size = 12)

print(p_line)
ggsave("plots/AtlasSD_vs_AtlasAlpha_LINEONLY.png", p_line, width = 7, height = 5, dpi = 300)

# optional: save the overlap table
fwrite(ov, "plots/overlap_AtlasSD_vs_AtlasAlpha.csv")

















# Assumes you already have `ov` with columns: alpha, SD_Atlas

# LINE ONLY (0–1 axes)
p_line_01 <- ggplot(ov, aes(x = alpha, y = SD_Atlas)) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Atlas SD vs Atlas alpha (overlapping CpGs) ",
    subtitle = sprintf("n=%d | ", nrow(ov)),
    x = "alpha ",
    y = "Standard_Deviation"
  ) +
  theme_minimal(base_size = 12)

print(p_line_01)
ggsave("plots/AtlasSD_vs_AtlasAlpha_LINEONLY_0to1.png", p_line_01, width = 7, height = 5, dpi = 300)


# POINTS + LINE (0–1 axes)
p_points_01 <- ggplot(ov, aes(x = alpha, y = SD_Atlas)) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Atlas SD vs Atlas alpha (overlapping CpGs) ",
    subtitle = sprintf("n=%d | ", nrow(ov)),
    x = "alpha ",
    y = "Standard_Deviation "
  ) +
  theme_minimal(base_size = 12)

print(p_points_01)
ggsave("plots/AtlasSD_vs_AtlasAlpha_POINTS_0to1.png", p_points_01, width = 7, height = 5, dpi = 300)






