
# Package intsallation
#-------------------------------------------------------------------------------    

install.packages("magrittr")
install.packages("cli")
install.packages("glue")   # to fix the permission denied one, after restart
install.packages("openxlsx")
install.packages("remotes")
remotes::install_github("al2na/methylKit")


# libraries 
#-------------------------------------------------------------------------------     

library(ggplot2)
library(methylKit)
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(scales)
  
#-------------------------------------------------------------------------------     

# Loading Bisulfite sequencing (methylation reads)
#-------------------------------------------------------------------------------     

Mdata <- read.xlsx("C:/Users/hamda/Desktop/UGI_thingy/data/gkac503_supplemental_files/SUPPLEMENTARY_TABLES_resubmission.xlsx", sheet = 6)

file1 <- "C:/Users/hamda/Desktop/UGI_thingy/3rdyear/genome/quality/methylation_hg38/Zvej16_bisulfite.trim_bismark_bt2.deduplicated.bismark.cov.gz"

Zvej16 <- methRead(
  location   = file1,
  sample.id  = "Zvej16",
  assembly   = "hg38",
  pipeline   = "bismarkCoverage",
  mincov     = 0,
  treatment  = 0      
)

View(Zvej16)


file2 <- "C:/Users/hamda/Desktop/UGI_thingy/3rdyear/genome/quality/methylation_hg38/SP75_bisulfite.trim_bismark_bt2.deduplicated.bismark.cov.gz"

SP75 <- methRead(
  location   = file2,
  sample.id  = "SP75",
  assembly   = "hg38",
  pipeline   = "bismarkCoverage",
  mincov     = 0,
  treatment  = 0      
)

SP75$methy <- SP75$numCs / SP75$coverage
Zvej16$methy <- Zvej16$numCs / Zvej16$coverage

file1 <- "C:/Users/hamda/Desktop/UGI_thingy/3rdyear/genome/quality/methylation_hg38/Zvej16_bisulfite.trim_bismark_bt2.deduplicated.bismark.cov.gz"

write.csv(Zvej16, file = "Zvej16file.csv")
write.csv(SP75, file = "SP75file.csv")

# fully loading the genomes 
#-------------------------------------------------------------------------------     
# Folder with your coverage files
cov_dir <- "C:/Users/hamda/Desktop/UGI_thingy/3rdyear/genome/bis/Bisgen/newtrim/methylation_hg38"

# Get all .bismark.cov or .bismark.cov.gz files
cov_files <- list.files(
  path       = cov_dir,
  pattern    = "\\.bismark\\.cov(\\.gz)?$",
  full.names = TRUE
)

# Make clean object names from filenames
sample_ids <- basename(cov_files)
sample_ids <- sub("_bisilfite.*", "", sample_ids)  # e.g. "Vac_164_Molar_bisilfite..." -> "Vac_164_Molar"

# Loop over files and create ONE data frame per file in the Global Environment
for (i in seq_along(cov_files)) {
  message("Reading: ", cov_files[i], "  -->  data frame: ", sample_ids[i])
  
  df <- read.table(
    cov_files[i],
    header = FALSE
  )
  
  # V7 = total coverage (methylated + unmethylated reads)
  df$V7 <- df$V5 + df$V6
  
  # Create an object with that name in your workspace
  assign(sample_ids[i], df, envir = .GlobalEnv)
}





#-------------------------------------------------------------------------------     



# load all the data into a list and count CgP appearances 
#-------------------------------------------------------------------------------     


# Load all CpG coverage tables and compute CpG presence across samples

# 1. Folder with your coverage files
cov_dir <- "C:/Users/hamda/Desktop/UGI_thingy/3rdyear/genome/bis/Bisgen/newtrim/methylation_hg38"

# 2. Get all .bismark.cov or .bismark.cov.gz files
cov_files <- list.files(
  path       = cov_dir,
  pattern    = "\\.bismark\\.cov(\\.gz)?$",
  full.names = TRUE
)

if (length(cov_files) == 0) {
  stop("No .bismark.cov or .bismark.cov.gz files found in the folder.")
}

# 3. Make clean object names from filenames
sample_ids <- basename(cov_files)
sample_ids <- sub("_bisilfite.*", "", sample_ids)  # e.g. "Vac_164_Molar_bisilfite..." -> "Vac_164_Molar"

# 4. Loop: read each file
#    - create ONE data frame per file in the Global Environment
#    - also store them in a named list 'cov_list'
cov_list <- vector("list", length(cov_files))
names(cov_list) <- sample_ids

for (i in seq_along(cov_files)) {
  message("Reading: ", cov_files[i], "  -->  data frame: ", sample_ids[i])
  
  df <- read.table(
    cov_files[i],
    header = FALSE
  )
  
  # V7 = total coverage (methylated + unmethylated reads)
  df$V7 <- df$V5 + df$V6
  
  # Put as its own object in your workspace
  assign(sample_ids[i], df, envir = .GlobalEnv)
  
  # Also keep in the list
  cov_list[[sample_ids[i]]] <- df
}

# At this point:
# - You have objects 10658, 11112, Vac_164_Molar, etc.
# - You also have cov_list[["10658"]], cov_list[["Vac_164_Molar"]], etc.

# 5. Compute CpG presence across samples
#    Define CpG by (chr, start) = (V1, V2)

cpg_counts_table <- table(unlist(
  lapply(cov_list, function(df) {
    # unique CpGs per sample to avoid double-counting within a sample
    unique(paste(df$V1, df$V2, sep = ":"))
  })
))

# Convert to data.frame
cpg_presence <- data.frame(
  CpG       = names(cpg_counts_table),
  n_samples = as.integer(cpg_counts_table),
  row.names = NULL
)

# Split CpG into chr and position
tmp <- do.call(rbind, strsplit(as.character(cpg_presence$CpG), ":", fixed = TRUE))
cpg_presence$chr <- tmp[, 1]
cpg_presence$pos <- as.integer(tmp[, 2])

# Reorder columns: chr, pos, n_samples
cpg_presence <- cpg_presence[, c("chr", "pos", "n_samples")]

# 6. Quick checks
# How many CpGs are seen in X samples?
table(cpg_presence$n_samples)

# First few CpGs with their presence counts
head(cpg_presence)

CpGpres <- as.data.frame(table(cpg_presence$n_samples))
#-------------------------------------------------------------------------------     




#-------------------------------------------------------------------------------     

# Exploratory analysis
#-------------------------------------------------------------------------------     

# density of methylation (for positions)

ggplot(Zvej16, aes(x=methy)) + 
  geom_density()

ggplot(SP75, aes(x=methy)) + 
  geom_density()


# histogram of methylation (for positions)
 
ggplot(Zvej16, aes(x=coverage)) + 
  geom_histogram()

ggplot(SP75, aes(x=methy)) + 
  geom_histogram()


# Presence of CpGs shared at 0/1 coverage 

ggplot( CpGpres , aes(x = Var1,  y = Freq)) + 
  geom()


# Forming a table which shows, sample, coverage, location etc
# may not be nessesary 
 

#-------------------------------------------------------------------------------     


# -------------------------------------------------------------------
# 0. Combine the list of coverage tables into one long data frame
#    Assumes cov_list already exists and each element has:
#    V1=chr, V2=pos, V5=methylated reads, V6=unmethylated reads, V7=coverage
# -------------------------------------------------------------------

# If V7 is not there for some reason, you can uncomment this block first:
# cov_list <- lapply(cov_list, function(df) {
#   if (!"V7" %in% names(df)) df$V7 <- df$V5 + df$V6
#   df
# })

trail_from_list <- imap_dfr(cov_list, ~{
  df <- .x
  tibble(
    chr   = df$V1,
    pos   = df$V2,
    ID    = .y,        # sample name from the list
    cov   = df$V7,     # total coverage at CpG
    mthyl = df$V5      # methylated reads at CpG
  )
})

# -------------------------------------------------------------------
# 1. All CpGs: "sites present at coverage ≥ X in at least N people"
# -------------------------------------------------------------------

cum_dist <- function(dat, min_cov) {
  site_counts <- dat %>%
    filter(!is.na(cov), cov >= min_cov) %>%           # cumulative: ≥ min_cov
    distinct(chr, pos, ID) %>%                        # one row per site/person
    count(chr, pos, name = "n_people")                # people per site
  
  dist <- site_counts %>%
    count(n_people, name = "n_sites") %>%             # sites with exactly N people
    arrange(n_people) %>%
    complete(n_people = full_seq(n_people, 1), fill = list(n_sites = 0)) %>%
    mutate(
      n_sites_at_least = rev(cumsum(rev(n_sites))),   # sites with ≥ N people
      n_people         = as.integer(n_people),
      min_cov          = min_cov,
      label            = paste0("≥", min_cov)
    )
  
  dist %>% select(n_people, n_sites_at_least, min_cov, label)
}

# Compute for thresholds 1..10 (min coverage back to 1)
dist_all_cum <- map_df(1:10, ~cum_dist(trail_from_list, .x)) %>%
  mutate(label = factor(label, levels = paste0("≥", 1:10)))

# Updated palette including ≥1
Hcolors <- c(
  "≥1"  = "grey80",
  "≥2"  = "grey60",
  "≥3"  = "grey40",
  "≥4"  = "grey20",
  "≥5"  = "#1b9e77",
  "≥6"  = "#d95f02",
  "≥7"  = "#7570b3",
  "≥8"  = "#e7298a",
  "≥9"  = "#66a61e",
  "≥10" = "#e6ab02"
)

pdf("C:/Users/hamda/Desktop/UGI_thingy/fig1_from_list_cov1.pdf", width = 10, height = 5)
ggplot() +
  geom_col(
    data = dist_all_cum,
    aes(x = n_people, y = n_sites_at_least, fill = label),
    position = "identity", width = 0.9, alpha = 0.65, colour = "black"
  ) +
  scale_fill_manual(values = Hcolors, name = "Min coverage") +
  scale_y_continuous(
    trans = "log10",
    labels = label_number(scale_cut = cut_short_scale()),
    breaks = 10^(0:6),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 12),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(
    title = "CpG sites shared by at least N people at different coverage thresholds",
    x = "N people (threshold)",
    y = "CpG sites (log scale)"
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        legend.position = "right")
dev.off()

# -------------------------------------------------------------------
# 2. Only methylated CpGs: "sites with methylation signal"
# -------------------------------------------------------------------

cum_dist_methylated <- function(dat, min_cov) {
  site_counts <- dat %>%
    filter(!is.na(cov), cov >= min_cov,
           !is.na(mthyl), mthyl > 0) %>%         # keep only methylated CpGs
    distinct(chr, pos, ID) %>%                   # one row per site/person
    count(chr, pos, name = "n_people")           # how many people have this site methylated
  
  dist <- site_counts %>%
    count(n_people, name = "n_sites") %>%        # sites with exactly N people methylated
    arrange(n_people) %>%
    complete(n_people = full_seq(n_people, 1), fill = list(n_sites = 0)) %>%
    mutate(
      n_sites_at_least = rev(cumsum(rev(n_sites))),  # sites with ≥ N people methylated
      n_people         = as.integer(n_people),
      min_cov          = min_cov,
      label            = paste0("≥", min_cov)
    )
  
  dist %>% select(n_people, n_sites_at_least, min_cov, label)
}

dist_all_cum_meth <- map_df(1:10, ~cum_dist_methylated(trail_from_list, .x)) %>%
  mutate(label = factor(label, levels = paste0("≥", 1:10)))

pdf("C:/Users/hamda/Desktop/UGI_thingy/fig2_from_list_cov1.pdf", width = 10, height = 5)
ggplot() +
  geom_col(
    data = dist_all_cum_meth,
    aes(x = n_people, y = n_sites_at_least, fill = label),
    position = "identity", width = 0.9, alpha = 0.65, colour = "black"
  ) +
  scale_fill_manual(values = Hcolors, name = "Min coverage") +
  scale_y_continuous(
    trans = "log10",
    labels = label_number(scale_cut = cut_short_scale()),
    breaks = 10^(0:6),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 12),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(
    title = "Methylated CpG sites shared by at least N people (by coverage threshold)",
    x = "N people",
    y = "Methylated CpG sites (log scale)"
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title.position = "plot",
        legend.position = "right")
dev.off()




# -------------------------------------------------------------------




# -------------------------------------------------------------------
# Cross genome methylation 


library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)

## 1. Collect all Vac data frames from the environment ------------------------

# Get object names that start with "Vac_"
vac_names <- ls(pattern = "^Vac_")

# Pull them into a named list
vac_list <- mget(vac_names)

# Check what we have
vac_names
# e.g. "Vac_164_Molar" "Vac_164_Petrous" ...

## 2. Build a long data frame: one row per CpG per sample ---------------------

vac_long <- imap_dfr(vac_list, ~{
  df <- .x
  
  tibble(
    sample    = .y,
    # bone type from the name
    bone      = case_when(
      str_detect(.y, "Molar")   ~ "Molar",
      str_detect(.y, "Petrous") ~ "Petrous",
      TRUE                      ~ "Other"
    ),
    # methylation percentage (column V4)
    perc_meth = df$V4
    # optional: filter out NAs if present
    # perc_meth = df$V4[!is.na(df$V4)]
  )
})

# Derive pair ID: "Vac_164_Molar" -> "Vac_164"
vac_long <- vac_long %>%
  mutate(pair = str_replace(sample, "_Molar|_Petrous", ""))

# Optional quick check:
# table(vac_long$pair, vac_long$bone)

## 3. Combined faceted density plot -------------------------------------------

# Choose colours for bone type
bone_cols <- c(Molar = "#1f77b4", Petrous = "#d62728", Other = "grey60")

p_all <- ggplot(vac_long, aes(x = perc_meth, colour = bone, fill = bone)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~ pair, scales = "free_y") +  # one panel per Vac pair
  scale_x_continuous(limits = c(0, 100)) +
  scale_colour_manual(values = bone_cols, name = "Bone") +
  scale_fill_manual(values = bone_cols, name = "Bone") +
  labs(
    title = "Distribution of CpG methylation percentages in Vac molar vs petrous bones",
    x = "Methylation percentage",
    y = "Density"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title.position = "plot"
  )

# Show in RStudio
print(p_all)

# Save combined figure
ggsave(
  filename = "C:/Users/hamda/Desktop/UGI_thingy/Vac_methyl_density_all_pairs.pdf",
  plot     = p_all,
  width    = 10,
  height   = 6
)

## 4. Separate density plot per pair (one PDF each) ---------------------------

pairs <- unique(vac_long$pair)

for (pr in pairs) {
  p_pair <- ggplot(filter(vac_long, pair == pr),
                   aes(x = perc_meth, colour = bone, fill = bone)) +
    geom_density(alpha = 0.3) +
    scale_x_continuous(limits = c(0, 100)) +
    scale_colour_manual(values = bone_cols, name = "Bone") +
    scale_fill_manual(values = bone_cols, name = "Bone") +
    labs(
      title = paste0("Methylation % density – ", pr),
      x = "Methylation percentage",
      y = "Density"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title.position = "plot"
    )
  
  ggsave(
    filename = paste0("C:/Users/hamda/Desktop/UGI_thingy/", pr, "_methyl_density.pdf"),
    plot     = p_pair,
    width    = 6,
    height   = 4
  )
}
















 

## 1. Collect all Vac data frames from the environment ------------------------

vac_names <- ls(pattern = "^Vac_")
vac_list  <- mget(vac_names)

# If V7 not present for some reason, uncomment this:
# vac_list <- lapply(vac_list, function(df) {
#   if (!"V7" %in% names(df)) df$V7 <- df$V5 + df$V6
#   df
# })

## 2. Build a long data frame: one row per CpG per sample ---------------------

vac_long <- imap_dfr(vac_list, ~{
  df <- .x
  tibble(
    sample    = .y,
    bone      = case_when(
      str_detect(.y, "Molar")   ~ "Molar",
      str_detect(.y, "Petrous") ~ "Petrous",
      TRUE                      ~ "Other"
    ),
    pair      = str_replace(.y, "_Molar|_Petrous", ""),
    cov       = df$V7,   # total coverage
    perc_meth = df$V4    # methylation percentage
  )
})

# Colours for bone type
bone_cols <- c(Molar = "#1f77b4", Petrous = "#d62728", Other = "grey60")

## Helper to make + save a faceted density plot ------------------------------

make_density_plot <- function(dat, title_suffix, filename_suffix) {
  p <- ggplot(dat, aes(x = perc_meth, colour = bone, fill = bone)) +
    geom_density(alpha = 0.3) +
    facet_wrap(~ pair, scales = "free_y") +
    scale_x_continuous(limits = c(0, 100)) +
    scale_colour_manual(values = bone_cols, name = "Bone") +
    scale_fill_manual(values = bone_cols, name = "Bone") +
    labs(
      title = paste0("Methylation % density – Vac molar vs petrous (", title_suffix, ")"),
      x = "Methylation percentage",
      y = "Density"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title.position = "plot"
    )
  
  print(p)
  
  ggsave(
    filename = paste0("C:/Users/hamda/Desktop/UGI_thingy/Vac_methyl_density_", filename_suffix, ".pdf"),
    plot     = p,
    width    = 10,
    height   = 6
  )
}

## 3. Plots for different minimum coverages ----------------------------------

# a) No coverage filter (all CpGs)
make_density_plot(
  dat            = vac_long,
  title_suffix   = "no coverage filter",
  filename_suffix = "cov_min0"
)

# b) Coverage ≥ 3
vac_long_cov3 <- vac_long %>% filter(cov >= 3)

make_density_plot(
  dat            = vac_long_cov3,
  title_suffix   = "coverage ≥ 3",
  filename_suffix = "cov_min3"
)

# c) Coverage ≥ 4
vac_long_cov4 <- vac_long %>% filter(cov >= 4)

make_density_plot(
  dat            = vac_long_cov4,
  title_suffix   = "coverage ≥ 4",
  filename_suffix = "cov_min4"
)





# --------------------------------------------------------------------------
# check between mdata and trail
#> trail_petrous$chrpos <- paste0(trail_petrous$chr, ":", trail_petrous$pos)

Mdata <- read.xlsx("data/gkac503_supplemental_files/SUPPLEMENTARY_TABLES_resubmission.xlsx", sheet = 6)

library(dplyr)

## 1. Prepare Mdata with hg38 key ----

Mdata_hg38 <- Mdata %>%
  # keep only rows that actually have hg38 coordinates
  filter(!is.na(chr_hg38), !is.na(hg38.pos)) %>%
  mutate(chrpos = paste0(chr_hg38, ":", hg38.pos))

## 2. Prepare petrous with the same key ----

trail_petrous2 <- trail_petrous %>%
  mutate(chrpos = paste0(chr, ":", pos))

## 3. CpGs shared between Mdata and petrous ----

shared_petrous <- Mdata_hg38 %>%
  inner_join(trail_petrous2, by = "chrpos")

# If you just want the list of unique CpG IDs:
shared_petrous_CpGs <- shared_petrous %>%
  distinct(CpG)

# Optional: peek
head(shared_petrous)
head(shared_petrous_CpGs)

## 4. Prepare molar with the same key ----
## (assuming you have a data frame called `trail_molar` with chr/pos/ID/cov/mthyl)

trail_molar2 <- trail_molar %>%
  mutate(chrpos = paste0(chr, ":", pos))

## 5. CpGs shared between Mdata and molar ----

shared_molar <- Mdata_hg38 %>%
  inner_join(trail_molar2, by = "chrpos")

shared_molar_CpGs <- shared_molar %>%
  distinct(CpG)

# Optional: peek
head(shared_molar)
head(shared_molar_CpGs)

