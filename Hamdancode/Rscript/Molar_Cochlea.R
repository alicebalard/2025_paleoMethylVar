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

# Example:
corr_for_pair("Vac_164", min_cov = 3)
corr_for_pair("Vac_210", min_cov = 10)

# -------------------------------------------------
 # 2A. Helper: turn a cov table into a vector of CpG keys #----------------------------------------------------------- get_cpg_keys <- function(df) { unique(paste(df$V1, df$V2, sep = ":")) } #----------------------------------------------------------- # 2B. Identify Vac pairs (Vac_164, Vac_XXX, ...) #----------------------------------------------------------- all_names <- names(cov_list) vac_names <- all_names[str_detect(all_names, "^Vac_")] vac_pairs <- vac_names %>% str_replace("_Molar|_Petrous", "") %>% unique() vac_pairs #----------------------------------------------------------- # 2C. Build overlap summary per pair #----------------------------------------------------------- vac_overlap_list <- map(vac_pairs, function(pr) { mol_name <- paste0(pr, "_Molar") pet_name <- paste0(pr, "_Petrous") # Skip if one member is missing if (!mol_name %in% all_names || !pet_name %in% all_names) return(NULL) mol_df <- cov_list[[mol_name]] pet_df <- cov_list[[pet_name]] mol_cpg <- get_cpg_keys(mol_df) pet_cpg <- get_cpg_keys(pet_df) shared <- intersect(mol_cpg, pet_cpg) mol_only <- setdiff(mol_cpg, pet_cpg) pet_only <- setdiff(pet_cpg, mol_cpg) tibble( pair = pr, n_molar = length(mol_cpg), n_petrous = length(pet_cpg), n_shared = length(shared), n_molar_only = length(mol_only), n_petrous_only = length(pet_only) ) }) vac_overlap_summary <- bind_rows(vac_overlap_list) vac_overlap_summary # Save summary table if you like write.csv( vac_overlap_summary, "C:/Users/hamda/Desktop/UGI_thingy/Vac_molar_petrous_CpG_overlap_summary.csv", row.names = FALSE )
