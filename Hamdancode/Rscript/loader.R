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
