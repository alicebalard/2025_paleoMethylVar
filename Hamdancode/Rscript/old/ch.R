# ------------------------------------------------------------------------------
#  chromomap_heatmap.R
#  Stand‑alone: draw an interactive CpG‑variance heat‑map with chromoMap
# ------------------------------------------------------------------------------
#  EXPECTED INPUT (already in your workspace)
#      var_cpg  : data.frame with columns
#                   Chr       integer 1‑22
#                   Pos       bp coordinate (hg38)
#                   var_beta  numeric value to visualise
#
#  USAGE
#      source("chromomap_heatmap.R")   # launches the map immediately
# ------------------------------------------------------------------------------

# === 1. PACKAGES =============================================================
# chromoMap uses base R colours unless you pass an explicit vector.  We rely on
# RColorBrewer for a clean sequential palette.
if (!requireNamespace("chromoMap",      quietly = TRUE)) install.packages("chromoMap")
if (!requireNamespace("readr",          quietly = TRUE)) install.packages("readr")
if (!requireNamespace("RColorBrewer",   quietly = TRUE)) install.packages("RColorBrewer")

library(chromoMap)
library(dplyr)
library(readr)
library(RColorBrewer)

# === 2. CHROMOSOME SIZES (hg38 autosomes) ====================================
hg38_sizes <- c(248956422, 242193529, 198295559, 190214555, 181538259,
                170805979, 159345973, 145138636, 138394717, 133797422,
                135086622, 133275309, 114364328, 107043718, 101991189,
                90338345,  83257441,  80373285,  58617616,  64444167,
                46709983,  50818468)

chrom_sizes <- data.frame(
  chrom = paste0("chr", 1:22),
  start = 1,
  end   = hg38_sizes
)

# === 3. PREPARE HEAT‑MAP DATA ===============================================
if (!exists("var_cpg")) stop("var_cpg not found – please load it before sourcing.")

hm <- var_cpg %>%
  mutate(id    = paste0("CpG_", row_number()),   # unique identifier for row‑names
         chrom = paste0("chr", Chr),
         start = Pos,
         end   = Pos,
         value = var_beta) %>%
  select(id, chrom, start, end, value)

# === 4. CREATE TEMPORARY TSVs FOR chromoMap ==================================
coords_file <- tempfile(fileext = ".txt")
data_file   <- tempfile(fileext = ".txt")
write_tsv(chrom_sizes, coords_file, col_names = FALSE)
write_tsv(hm,           data_file,   col_names = FALSE)

# === 5. PLOT IMMEDIATELY =====================================================
# Use a Brewer sequential palette (9 shades of blue) — convert to hex vector
blue_ramp <- brewer.pal(9, "Blues")

chromoMap(
  coords_file,
  data_file,
  data_type     = "numeric",            # heat‑map
  data_colors   = list(blue_ramp),       # vector of hex colours
  legend        = TRUE,
  canvas_width  = 900,
  canvas_height = 500
)

# ------------------------------------------------------------------------------
#  END OF FILE
# ------------------------------------------------------------------------------
