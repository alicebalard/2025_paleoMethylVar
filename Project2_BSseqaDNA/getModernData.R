## Script to get the full modern data results with pr(hv) calculated
## Alice Balard, 2 dec 2025
###########################
library(progress)
library(data.table)

## On bash, clone my repo
## git clone git@github.com:alicebalard/2024_hvCpG.git

# Define parent folder containing all "Atlas_batchXXX" folders (put the correct full path)
parent_dir <- "~/Documents/GIT/2024_hvCpG/05_hvCpGalgorithm/resultsDir/Atlas/Atlas10X/"

prepAtlasdt <- function(){

  # Get list of relevant RDS files
  rds_files <- dir(parent_dir, pattern = "results_Atlas10X_[0-9]+CpGs_0_8p0_0_65p1\\.rds$", 
                   recursive = TRUE, full.names = TRUE)
  
  all_cpg_values <- numeric()
  pb <- progress_bar$new(total = length(rds_files), format = "📦 :current/:total [:bar] :percent")
  
  for (file in rds_files) {
    obj <- readRDS(file)
    if (is.matrix(obj)) obj <- obj[, 1]
    all_cpg_values <- c(all_cpg_values, obj)
    pb$tick()
  }
  
  # Create data.table from named vector
  dt <- data.table(
    name =   sub("-[0-9]+$", "", names(all_cpg_values)), # just keep the C position instead of C + 1
    alpha = as.numeric(all_cpg_values)
  )
  
  #######################################################################
  # Parse "chr_pos" in name into chr, start_pos. NB: takes a couple of minutes
  dt[, c("chr", "pos") := tstrsplit(name, "_", fixed = TRUE)]
  
  # Convert to integer/numeric if not already
  dt[, pos := as.integer(pos)]
  
  ## Check chromosomes present:
  message("Chromosomes in the dataset:")
  table(unique(dt$chr))
  
  # Convert chr from "chrN" to  factor
  dt[, chr := sub("chr", "", chr)]
  dt[, chr := factor(chr, levels = as.character(c(1:22, "X", "Y", "M")))]
  
  ## Check chromosomes order:
  message("Chromosomes in the dataset:")
  table(unique(dt$chr))
  
  # Compute cumulative position offsets for Manhattan plot
  setorder(dt, chr, pos)
  
  offsets <- dt[, .(max_pos = max(pos, na.rm = TRUE)), by = chr]
  offsets[, cum_offset := c(0, head(cumsum(as.numeric(max_pos)), -1))]
  
  dt <- merge(dt, offsets[, .(chr, cum_offset)], by = "chr", all.x = TRUE, sort = FALSE)
  
  # Convert to integer/numeric if not already
  dt[, cum_offset := as.numeric(cum_offset)]
  dt[, pos2 := pos + cum_offset]
  
  return(dt)
}

system.time(Atlas_dt <- prepAtlasdt())

nrow(Atlas_dt) # 23036026
head(Atlas_dt)
