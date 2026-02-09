# ------------------------------------------------------------------------------
#  methylation_variance_analysis.R
#  End‑to‑end pipeline for CpG‑methylation variance + visualisation
#  Author: <your‑name>   Last update: 2025‑06‑17
# ------------------------------------------------------------------------------

# === 0. SETUP ================================================================
required <- c("dplyr", "ggplot2", "tidyr", "scales", "forcats", "readr", "chromoMap")
new_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(new_pkgs)) install.packages(new_pkgs)
invisible(lapply(required, library, character.only = TRUE))

# -----------------------------------------------------------------------------
# 1. TIDY RAW DATA -------------------------------------------------------------
#    Assumes a raw data.frame `testmerge` with columns:
#      * V1            = "chr_position" (e.g. "4_91523712")
#      * V5            = β‑value
#      * sample_alias  = individual ID
#      * period4ana    = archaeological period
#    The helper returns a clean tibble with:
#      Chr (int), Pos (int), beta (dbl), sample_alias, period4ana
# -----------------------------------------------------------------------------

tidy_cpg <- function(df, keep_autosomes = TRUE) {
  df %>%
    tidyr::separate(V1, into = c("Chr", "Pos"), sep = "_", convert = TRUE) %>%
    dplyr::distinct(sample_alias, Chr, Pos, period4ana, .keep_all = TRUE) %>%
    dplyr::filter(!is.na(V5)) %>%
    {
      if (keep_autosomes) dplyr::filter(., Chr %in% 1:22) else .
    } %>%
    dplyr::rename(beta = V5) %>%
    dplyr::arrange(Chr, Pos)
}

# Cleaned data ---------------------------------------------------------------
tm2 <- tidy_cpg(testmerge)

# -----------------------------------------------------------------------------
# 2. VARIANCE PER CHR × PERIOD  (CpGs as replicates) ---------------------------
# -----------------------------------------------------------------------------

var_chr_period <- tm2 %>%
  dplyr::group_by(period4ana, Chr) %>%
  dplyr::summarise(var_beta = stats::var(beta, na.rm = TRUE),
                   n_sites  = dplyr::n(),
                   .groups  = "drop")

# -----------------------------------------------------------------------------
# 3. BAR‑PLOT OF CHROMOSOMAL VARIANCE (FACETED BY PERIOD) ----------------------
# -----------------------------------------------------------------------------

plot_chr_variance <- function(data = var_chr_period) {
  ggplot2::ggplot(data,
                  ggplot2::aes(x = forcats::fct_inorder(factor(Chr)),
                               y = var_beta)) +
    ggplot2::geom_col(fill = "#2ca25f") +
    ggplot2::facet_wrap(~ period4ana, nrow = 1) +
    ggplot2::labs(x = "Chromosome", y = "Variance of methylation (β)") +
    ggplot2::theme_minimal()
}

# -----------------------------------------------------------------------------
# 4. VARIANCE PER SITE  --------------------------------------------------------
# -----------------------------------------------------------------------------

site_var <- tm2 %>%
  dplyr::group_by(Chr, Pos, period4ana) %>%
  dplyr::filter(dplyr::n() > 1) %>%
  dplyr::summarise(var_site = stats::var(beta, na.rm = TRUE), .groups = "drop")

var_cpg <- site_var %>%
  dplyr::group_by(Chr, Pos) %>%
  dplyr::summarise(var_beta = mean(var_site), .groups = "drop")

# -----------------------------------------------------------------------------
# 5. SCATTER PLOT ALONG CHROMOSOMES  ------------------------------------------
# -----------------------------------------------------------------------------

plot_scatter <- function(data = var_cpg, top_prop = 0.01) {
  q_cut <- stats::quantile(data$var_beta, 1 - top_prop, na.rm = TRUE)
  ggplot2::ggplot(data,
                  ggplot2::aes(Pos / 1e6, var_beta,
                               colour = var_beta > q_cut)) +
    ggplot2::geom_point(alpha = 0.6, size = 0.7) +
    ggplot2::scale_colour_manual(values = c("grey60", "red"), guide = "none") +
    ggplot2::facet_wrap(~ Chr, scales = "free_x", ncol = 4) +
    ggplot2::scale_x_continuous(labels = scales::label_number(suffix = " Mb")) +
    ggplot2::labs(x = "Genomic position (Mb)",
                  y = "Variance of methylation (β)") +
    ggplot2::theme_bw()
}

# -----------------------------------------------------------------------------
# 6. INTERACTIVE HEATMAP (chromoMap) ------------------------------------------
# -----------------------------------------------------------------------------

plot_chromomap <- function(data = var_cpg) {
  hg38_sizes <- c(248956422, 242193529, 198295559, 190214555, 181538259,
                  170805979, 159345973, 145138636, 138394717, 133797422,
                  135086622, 133275309, 114364328, 107043718, 101991189,
                  90338345, 83257441, 80373285, 58617616, 64444167,
                  46709983, 50818468)
  
  chrom_sizes <- data.frame(
    chrom = paste0("chr", 1:22),
    start = 1,
    end   = hg38_sizes
  )
  
  hm <- data %>%
    dplyr::mutate(chrom = paste0("chr", Chr),
                  start = Pos,
                  end   = Pos,
                  value = var_beta) %>%
    dplyr::select(chrom, start, end, value)
  
  coords_file <- tempfile(fileext = ".txt")
  data_file   <- tempfile(fileext = ".txt")
  
  readr::write_tsv(chrom_sizes, coords_file, col_names = FALSE)
  readr::write_tsv(hm,           data_file,   col_names = FALSE)
  
  chromoMap::chromoMap(
    coords_file,
    data_file,
    data_type     = "heatmap",
    legend        = TRUE,
    cmap          = "Blues",
    aggr_type     = "average",
    canvas_width  = 900,
    canvas_height = 500,
    exportoptions = TRUE
  )
}

# -----------------------------------------------------------------------------
# 7. SHOW / SAVE PLOTS ---------------------------------------------------------
# -----------------------------------------------------------------------------

# Helper to both display on‑screen *and* write a PDF if desired
save_and_show <- function(plot_expr, filename = NULL, width = 8, height = 6) {
  p <- rlang::eval_tidy(rlang::enquo(plot_expr))
  print(p)                                   # show immediately
  if (!is.null(filename)) {
    ggplot2::ggsave(filename, plot = p, width = width, height = height)
  }
  invisible(p)
}

if (interactive()) {
  save_and_show(plot_chr_variance(), "chr_variance_bar.pdf", width = 8, height = 3)
  save_and_show(plot_scatter(),       "cpg_scatter_facets.pdf", width = 9, height = 7)
  plot_chromomap()  # opens interactive viewer; use toolbar to export PNG/SVG
}


if (interactive()) {
  if ("bar"       %in% RUN_PLOTS) save_and_show(plot_chr_variance(),
                                                "chr_variance_bar.pdf",
                                                width = 8, height = 3)
  if ("scatter"   %in% RUN_PLOTS) save_and_show(plot_scatter(),
                                                "cpg_scatter_facets.pdf",
                                                width = 9, height = 7)
  if ("chromomap" %in% RUN_PLOTS) plot_chromomap()
}


# === END =====================================================================
