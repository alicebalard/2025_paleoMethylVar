library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(tibble)

#-----------------------------------------------------------
# 1A. Identify Vac molar / Vac petrous / other samples
#-----------------------------------------------------------

all_names  <- names(cov_list)
vac_names  <- all_names[str_detect(all_names, "^Vac_")]

vac_molar   <- vac_names[str_detect(vac_names, "Molar")]
vac_petrous <- vac_names[str_detect(vac_names, "Petrous")]

other_samples <- setdiff(all_names, vac_names)

vac_molar
vac_petrous
other_samples

#-----------------------------------------------------------
# 1B. Build two variant lists
#     - cov_list_molar   = all non-Vac + Vac_*_Molar only
#     - cov_list_petrous = all non-Vac + Vac_*_Petrous only
#-----------------------------------------------------------

cov_list_molar   <- cov_list[c(vac_molar)] # made so that it will only take the molar samples 
cov_list_petrous <- cov_list[c(other_samples, vac_petrous)]
 
#-----------------------------------------------------------
# 1C. Helper to convert a list of coverage tables into one long df
#-----------------------------------------------------------

make_trail <- function(clist) {
  purrr::imap_dfr(clist, ~{
    df <- .x
    tibble(
      chr   = df$V1,
      pos   = df$V2,
      ID    = .y,
      cov   = df$V7,  # total coverage
      mthyl = df$V5   # methylated reads
    )
  })
}

# "All samples" version (keeps both bones for Vac)
trail_all     <- make_trail(cov_list)

# "Vac molar only" + all others
trail_molar   <- make_trail(cov_list_molar)

# "Vac petrous only" + all others
trail_petrous <- make_trail(cov_list_petrous)

# Quick checks
head(trail_molar)
head(trail_petrous)
table(trail_molar$ID %>% str_detect("^Vac_"))
table(trail_petrous$ID %>% str_detect("^Vac_"))
