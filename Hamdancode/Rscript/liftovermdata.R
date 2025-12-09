library(dplyr)

# Make sure hg19.pos is numeric
Mdata$hg19.pos <- as.integer(as.character(Mdata$hg19.pos))

# BED is 0-based, half-open: start = pos-1, end = pos
Mdata_bed <- Mdata %>%
  transmute(
    chr  = chr,                 # e.g. "chr1"
    start = hg19.pos - 1L,      # 0-based
    end   = hg19.pos,           # 1-based end
    name  = CpG                 # CpG ID as name
  )

# Write to your Desktop folder
bed_path <- "C:/Users/hamda/Desktop/UGI_thingy/Mdata_hg19.bed"
write.table(
  Mdata_bed,
  file = bed_path,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
bed_path


library(dplyr)

# Read lifted hg38 BED
lifted <- read.table(
  "C:/Users/hamda/Desktop/UGI_thingy/Mdata_hg38.bed",
  sep = "\t",
  stringsAsFactors = FALSE,
  col.names = c("chr_hg38", "start_hg38", "end_hg38", "CpG")
)

# Merge hg38 coords back into Mdata by CpG ID
Mdata <- Mdata %>%
  left_join(
    lifted %>%
      transmute(
        CpG,
        chr_hg38,
        hg38.pos = end_hg38   # 1-based position
      ),
    by = "CpG"
  )

head(Mdata[, c("CpG", "chr", "hg19.pos", "chr_hg38", "hg38.pos")])

Mdata$chr_hg38 <- Mdata$chr_hg38.y
Mdata$hg38.pos <- Mdata$hg38.pos.y

library(dplyr)

## 1. Prepare Mdata hg38 key ----

Mdata_hg38 <- Mdata %>%
  filter(!is.na(chr_hg38), !is.na(hg38.pos)) %>%    # keep only mapped CpGs
  mutate(chrpos = paste0(chr_hg38, ":", hg38.pos))

## 2. Prepare petrous with the same key ----

trail_petrous2 <- trail_petrous %>%
  mutate(chrpos = paste0(chr, ":", pos))

## 3. CpGs present in both Mdata and petrous ----

shared_petrous <- Mdata_hg38 %>%
  inner_join(trail_petrous2, by = "chrpos")

# Unique CpG IDs that are shared:
shared_petrous_CpGs <- shared_petrous %>%
  distinct(CpG)

## 4. Prepare molar with the same key ----

trail_molar2 <- trail_molar %>%
  mutate(chrpos = paste0(chr, ":", pos))

## 5. CpGs present in both Mdata and molar ----

shared_molar <- Mdata_hg38 %>%
  inner_join(trail_molar2, by = "chrpos")

shared_molar_CpGs <- shared_molar %>%
  distinct(CpG)

## 6. CpGs present in ALL THREE: Mdata, petrous, molar (optional) ----

CpGs_in_all_three <- intersect(shared_petrous_CpGs$CpG,
                               shared_molar_CpGs$CpG)

# Inspect:
head(shared_petrous)
head(shared_molar)
head(shared_petrous_CpGs)
head(shared_molar_CpGs)
CpGs_in_all_three



 