library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(tibble)

# 1. Folder with your coverage files
cov_dir <- "C:/Users/hamda/Desktop/UGI_thingy/3rdyear/genome/bis/Bisgen/newtrim/methylation_hg38"

# 2. List all .bismark.cov and .bismark.cov.gz files
cov_files <- list.files(
  path       = cov_dir,
  pattern    = "\\.bismark\\.cov(\\.gz)?$",
  full.names = TRUE
)

# 3. Clean sample IDs (handles 'bisilfite' and 'bisulfite')
sample_ids <- basename(cov_files)
sample_ids <- sub("_bis[ui]lfite.*", "", sample_ids)



# 4. Read each file, add V7, assign + store in list
cov_list <- vector("list", length(cov_files))
names(cov_list) <- sample_ids

for (i in seq_along(cov_files)) {
  message("Reading: ", cov_files[i], "  -->  data frame: ", sample_ids[i])
  
  df <- read.table(cov_files[i], header = FALSE)
  df$V7 <- df$V5 + df$V6   # total coverage
  
  assign(sample_ids[i], df, envir = .GlobalEnv)
  cov_list[[sample_ids[i]]] <- df
}



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

 



 