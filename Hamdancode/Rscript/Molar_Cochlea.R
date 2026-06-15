library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(tidyr)

library(dplyr)
library(purrr)
library(stringr)
library(tibble)

#-----------------------------------------------------------
# 2A. Helper: turn a cov table into a vector of CpG keys
#-----------------------------------------------------------

get_cpg_keys <- function(df) {
  unique(paste(df$V1, df$V2, sep = ":"))
}

#-----------------------------------------------------------
# 2B. Identify Vac pairs (Vac_164, Vac_XXX, ...)
#-----------------------------------------------------------

all_names <- names(cov_list)
vac_names <- all_names[str_detect(all_names, "^Vac_")]

vac_pairs <- vac_names %>%
  str_replace("_Molar|_Petrous", "") %>%
  unique()

vac_pairs

#-----------------------------------------------------------
# 2C. Build overlap summary per pair
#-----------------------------------------------------------

vac_overlap_list <- map(vac_pairs, function(pr) {
  mol_name  <- paste0(pr, "_Molar")
  pet_name  <- paste0(pr, "_Petrous")
  
  # Skip if one member is missing
  if (!mol_name %in% all_names || !pet_name %in% all_names) return(NULL)
  
  mol_df <- cov_list[[mol_name]]
  pet_df <- cov_list[[pet_name]]
  
  mol_cpg <- get_cpg_keys(mol_df)
  pet_cpg <- get_cpg_keys(pet_df)
  
  shared      <- intersect(mol_cpg, pet_cpg)
  mol_only    <- setdiff(mol_cpg, pet_cpg)
  pet_only    <- setdiff(pet_cpg, mol_cpg)
  
  tibble(
    pair            = pr,
    n_molar         = length(mol_cpg),
    n_petrous       = length(pet_cpg),
    n_shared        = length(shared),
    n_molar_only    = length(mol_only),
    n_petrous_only  = length(pet_only)
  )
})

vac_overlap_summary <- bind_rows(vac_overlap_list)
vac_overlap_summary

# Save summary table if you like
write.csv(
  vac_overlap_summary,
  "C:/Users/hamda/Desktop/UGI_thingy/Vac_molar_petrous_CpG_overlap_summary.csv",
  row.names = FALSE
)


#-----------------------------------------------------------
# 2D. CpG lists per pair, with methylation % AND coverage
#     vac_overlap_tables[["Vac_164"]]$shared has:
#     chr, pos, meth_molar, cov_molar, meth_petrous, cov_petrous
#-----------------------------------------------------------

vac_overlap_tables <- map(vac_pairs, function(pr) {
  mol_name  <- paste0(pr, "_Molar")
  pet_name  <- paste0(pr, "_Petrous")
  
  if (!mol_name %in% all_names || !pet_name %in% all_names) return(NULL)
  
  mol_df <- cov_list[[mol_name]]
  pet_df <- cov_list[[pet_name]]
  
  # If coverage column V7 somehow missing, rebuild it
  if (!"V7" %in% names(mol_df)) mol_df$V7 <- mol_df$V5 + mol_df$V6
  if (!"V7" %in% names(pet_df)) pet_df$V7 <- pet_df$V5 + pet_df$V6
  
  # CpG keys ("chr:pos")
  mol_cpg <- get_cpg_keys(mol_df)
  pet_cpg <- get_cpg_keys(pet_df)
  
  shared      <- intersect(mol_cpg, pet_cpg)
  mol_only    <- setdiff(mol_cpg, pet_cpg)
  pet_only    <- setdiff(pet_cpg, mol_cpg)
  
  # Helper: turn "chr:pos" vector into df
  split_to_df <- function(vec) {
    if (length(vec) == 0) return(tibble(chr = character(), pos = integer()))
    tmp <- str_split_fixed(vec, ":", 2)
    tibble(
      chr = tmp[, 1],
      pos = as.integer(tmp[, 2])
    )
  }
  
  # Collapse methylation + coverage per CpG in each sample
  # V4 = % methylation, V7 = coverage
  mol_info <- mol_df %>%
    transmute(chr = V1, pos = V2,
              meth_molar = V4,
              cov_molar  = V7) %>%
    group_by(chr, pos) %>%
    summarise(
      meth_molar = mean(meth_molar, na.rm = TRUE),
      cov_molar  = mean(cov_molar,  na.rm = TRUE),
      .groups    = "drop"
    )
  
  pet_info <- pet_df %>%
    transmute(chr = V1, pos = V2,
              meth_petrous = V4,
              cov_petrous  = V7) %>%
    group_by(chr, pos) %>%
    summarise(
      meth_petrous = mean(meth_petrous, na.rm = TRUE),
      cov_petrous  = mean(cov_petrous,  na.rm = TRUE),
      .groups      = "drop"
    )
  
  # Shared CpGs with methylation % and coverage from both bones
  shared_df <- split_to_df(shared) %>%
    left_join(mol_info, by = c("chr", "pos")) %>%
    left_join(pet_info, by = c("chr", "pos"))
  
  # Molar-only CpGs (with their methylation + coverage in molar)
  molar_only_df <- split_to_df(mol_only) %>%
    left_join(mol_info, by = c("chr", "pos"))
  
  # Petrous-only CpGs (with their methylation + coverage in petrous)
  petrous_only_df <- split_to_df(pet_only) %>%
    left_join(pet_info, by = c("chr", "pos"))
  
  list(
    shared       = shared_df,       # chr, pos, meth_molar, cov_molar, meth_petrous, cov_petrous
    molar_only   = molar_only_df,   # chr, pos, meth_molar, cov_molar
    petrous_only = petrous_only_df  # chr, pos, meth_petrous, cov_petrous
  )
})

names(vac_overlap_tables) <- vac_pairs

View(vac_overlap_tables[["Vac_164"]]$shared)



#Correlation checking thing 
#------------------------------------------------------------------
# Checking for correlation between molar and petrous in Vac 164

corr_for_pair <- function(pair_id, min_cov = 3) {
  df <- vac_overlap_tables[[pair_id]]$shared %>%
    filter(
      !is.na(meth_molar),
      !is.na(meth_petrous),
      !is.na(cov_molar),
      !is.na(cov_petrous),
      cov_molar   >= min_cov,
      cov_petrous >= min_cov
    )
  
  message("Pair: ", pair_id,
          " | min_cov = ", min_cov,
          " | n CpGs = ", nrow(df))
  
  print(cor.test(df$meth_molar, df$meth_petrous, method = "pearson"))
  
  p <- ggplot(df, aes(meth_molar, meth_petrous)) +
    geom_point(alpha = 0.3, size = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_smooth(method = "lm", se = FALSE) +
    coord_fixed(xlim = c(0, 100), ylim = c(0, 100)) +
    labs(
      x = "Molar methylation (%)",
      y = "Petrous methylation (%)",
      title = paste0("CpG-wise methylation correlation (",
                     pair_id, ", min cov = ", min_cov, ")")
    ) +
    theme_minimal()
  
  print(p)
  invisible(list(data = df, plot = p))
}

 



suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(stringr)
  library(tibble)
})

# -----------------------------
# SETTINGS
# -----------------------------
min_cov <- 3
out_dir <- "plots/VAC_pair_correlations"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

safe_file <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("[^A-Za-z0-9._-]+", "_")
}

# -----------------------------
# Robust per-pair runner
# -----------------------------
run_corr_pair <- function(pair_id, min_cov = 3, out_dir = "plots/VAC_pair_correlations") {
  
  obj <- vac_overlap_tables[[pair_id]]
  
  # default summary row (filled in if we can compute correlation)
  out <- tibble(
    pair_id = pair_id,
    min_cov = min_cov,
    status  = NA_character_,
    n_cpgs  = NA_integer_,
    pearson_r = NA_real_,
    pearson_p = NA_real_,
    pearson_t = NA_real_,
    pearson_df = NA_real_,
    pearson_ci_low = NA_real_,
    pearson_ci_high = NA_real_,
    plot_path = NA_character_
  )
  
  plot_path_png <- file.path(out_dir, paste0("VAC_", safe_file(pair_id), "_mincov", min_cov, ".png"))
  out$plot_path <- plot_path_png
  
  # -------- handle NULL / missing shared ----------
  if (is.null(obj) || is.null(obj$shared)) {
    out$status <- "shared_is_NULL"
    out$n_cpgs <- 0L
    
    p_empty <- ggplot() +
      annotate("text", x = 0, y = 0, label = "NO DATA: shared table is NULL", size = 5) +
      labs(
        title = paste0("CpG-wise methylation correlation (", pair_id, ", min cov = ", min_cov, ")")
      ) +
      theme_void()
    
    ggsave(plot_path_png, p_empty, width = 7.5, height = 5.5, dpi = 300)
    return(list(summary = out, plot = p_empty, data = NULL))
  }
  
  df0 <- obj$shared
  
  # -------- handle weird types ----------
  if (!is.data.frame(df0)) {
    out$status <- paste0("shared_not_df: ", class(df0)[1])
    out$n_cpgs <- 0L
    
    p_empty <- ggplot() +
      annotate("text", x = 0, y = 0, label = "NO DATA: shared is not a data.frame", size = 5) +
      labs(
        title = paste0("CpG-wise methylation correlation (", pair_id, ", min cov = ", min_cov, ")")
      ) +
      theme_void()
    
    ggsave(plot_path_png, p_empty, width = 7.5, height = 5.5, dpi = 300)
    return(list(summary = out, plot = p_empty, data = NULL))
  }
  
  # -------- do the filtering ----------
  df <- df0 %>%
    filter(
      !is.na(meth_molar),
      !is.na(meth_petrous),
      !is.na(cov_molar),
      !is.na(cov_petrous),
      cov_molar   >= min_cov,
      cov_petrous >= min_cov
    )
  
  n_cpgs <- nrow(df)
  out$n_cpgs <- n_cpgs
  
  # -------- plot ----------
  if (n_cpgs == 0) {
    out$status <- "no_CpGs_after_filter"
    
    p_empty <- ggplot() +
      annotate("text", x = 0, y = 0,
               label = paste0("NO CpGs after filtering\n(min_cov=", min_cov, ")"),
               size = 5) +
      labs(
        title = paste0("CpG-wise methylation correlation (", pair_id, ", min cov = ", min_cov, ")")
      ) +
      theme_void()
    
    ggsave(plot_path_png, p_empty, width = 7.5, height = 5.5, dpi = 300)
    return(list(summary = out, plot = p_empty, data = df))
  }
  
  p <- ggplot(df, aes(meth_molar, meth_petrous)) +
    geom_point(alpha = 0.3, size = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_smooth(method = "lm", se = FALSE) +
    coord_fixed(xlim = c(0, 100), ylim = c(0, 100)) +
    labs(
      x = "Molar methylation (%)",
      y = "Petrous methylation (%)",
      title = paste0("CpG-wise methylation correlation (", pair_id, ", min cov = ", min_cov, ")"),
      subtitle = paste0("n CpGs = ", n_cpgs)
    ) +
    theme_minimal(base_size = 13)
  
  ggsave(plot_path_png, p, width = 7.5, height = 5.5, dpi = 300)
  
  # -------- correlation (guard against zero variance / too few points) ----------
  if (n_cpgs >= 3 &&
      sd(df$meth_molar, na.rm = TRUE) > 0 &&
      sd(df$meth_petrous, na.rm = TRUE) > 0) {
    
    ct <- tryCatch(
      cor.test(df$meth_molar, df$meth_petrous, method = "pearson"),
      error = function(e) NULL
    )
    
    if (!is.null(ct)) {
      out$status <- "ok"
      out$pearson_r <- unname(ct$estimate)
      out$pearson_p <- unname(ct$p.value)
      out$pearson_t <- unname(ct$statistic)
      out$pearson_df <- unname(ct$parameter)
      out$pearson_ci_low  <- unname(ct$conf.int[1])
      out$pearson_ci_high <- unname(ct$conf.int[2])
    } else {
      out$status <- "cor_test_failed"
    }
  } else {
    out$status <- "too_few_or_zero_variance"
  }
  
  list(summary = out, plot = p, data = df)
}

# -----------------------------
# Run ALL pairs (won't crash now)
# -----------------------------
pair_ids <- names(vac_overlap_tables)

results_list <- purrr::map(pair_ids, ~ run_corr_pair(.x, min_cov = min_cov, out_dir = out_dir))
names(results_list) <- pair_ids

# -----------------------------
# Summary table (1 row per pair)
# -----------------------------
corr_summary <- bind_rows(purrr::map(results_list, "summary")) %>%
  arrange(desc(n_cpgs))

print(corr_summary)

write.csv(
  corr_summary,
  file = file.path(out_dir, paste0("VAC_correlation_summary_mincov", min_cov, ".csv")),
  row.names = FALSE
)

cat("Saved", nrow(corr_summary), "plots to:", out_dir, "\n")
