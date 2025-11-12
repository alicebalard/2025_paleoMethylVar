## From: https://onlinelibrary.wiley.com/doi/10.1111/eva.13743
## "We called all CG dinucleotide autosomal positions (n=26,752,702) from the human (hg19) reference genome using the R Bioconductor package “BSgenome.Hsapiens.UCSC.hg19” (Pagès, 2019) and stored these in a BED file. We then filtered these by removing any positions overlapping with SNP positions from dbSNP 142 (Sherry et al., 2001). Our aim here was to avoid confounding between methylation signals and real variants at CpG positions. There remained 13,270,411 autosomal CpG positions in the reference genome."

## if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
## BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg19", "GenomicRanges", "rtracklayer",
##                       "SNPlocs.Hsapiens.dbSNP142.GRCh37"))

library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(Biostrings)
library(rtracklayer)
library(SNPlocs.Hsapiens.dbSNP142.GRCh37)
library(BiocParallel)

setwd("/SAN/ghlab/pophistory/Alice/paleo_project/data/03refBed")

# Step 1: Extract all autosomal CpG positions (BSgenome: "chr1", "chr2", ...)
autosomes <- paste0("chr", 1:22)
genome <- BSgenome.Hsapiens.UCSC.hg19

get_cpg_gr <- function(chr) {
  hits <- matchPattern("CG", genome[[chr]])
  if (length(hits) > 0) {
    GRanges(seqnames=chr, ranges=IRanges(start=start(hits), width=2))
  } else {
    GRanges()
  }
}

register(MulticoreParam(workers=4))
cpg_list <- bplapply(autosomes, get_cpg_gr)
cpg_gr <- do.call(c, cpg_list)
cat("Total CpGs expected: 26752702", "\n")
cat("Total CpGs found:", length(cpg_gr), "\n")
export(cpg_gr, "all_autosomal_CpGs_hg19.bed", format="BED")

# Step 2: Extract all dbSNP142 SNPs for autosomes (SNPlocs: "ch1", "ch2", ...)
snps <- SNPlocs.Hsapiens.dbSNP142.GRCh37
autosomes_snps <- paste0("ch", 1:22)

get_snp_gr <- function(chr) {
  snp_chr <- snpsBySeqname(snps, chr)
  if (length(snp_chr) > 0) {
    GRanges(seqnames=chr, ranges=IRanges(start=start(snp_chr), width=1))
  } else {
    GRanges()
  }
}

snp_list <- bplapply(autosomes_snps, get_snp_gr)
snp_gr <- do.call(c, snp_list)
cat("Total dbSNP142 SNPs:", length(snp_gr), "\n")

# Step 3: Remove CpG positions overlapping dbSNP142 SNPs
# (need to harmonize seqnames: change "chr1" to "ch1" in cpg_gr or vice versa)
seqlevels(cpg_gr) <- gsub("chr", "ch", seqlevels(cpg_gr))
seqnames(cpg_gr) <- gsub("chr", "ch", seqnames(cpg_gr))

overlaps <- findOverlaps(cpg_gr, snp_gr, ignore.strand=TRUE)
cpg_filtered <- cpg_gr[-queryHits(overlaps)]
cat("Filtered CpGs:", length(cpg_filtered), "\n")
export(cpg_filtered, "autosomal_CpGs_hg19_no_dbSNP142.bed", format="BED")

cat("CpG positions saved to 'hg19_all_cpg.bed'\n")

























## Install required Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg19", "GenomicRanges", "rtracklayer"))

library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(rtracklayer)

# Extract all chromosomes (including autosomes, X, Y, and mitochondrial DNA)
all_chr <- seqnames(BSgenome.Hsapiens.UCSC.hg19)

# Initialize an empty GRanges object for CpG positions
cpg_positions <- GRanges()

# Iterate over each chromosome to find CpG positions
for (chr in all_chr) {
    cat("Processing", chr, "\n")
    # Get the sequence for the chromosome
    chr_seq <- BSgenome.Hsapiens.UCSC.hg19[[chr]]
    
    # Find CpG dinucleotide positions
    cpg_matches <- matchPattern("CG", chr_seq)
    
    # Store positions as GRanges object
    cpg_positions <- c(cpg_positions, GRanges(
        seqnames = chr,
        ranges = IRanges(start = start(cpg_matches), width = 2),
        strand = "*"
    ))
}

# Sort and remove duplicates
cpg_positions <- sort(unique(cpg_positions))

# Export CpG positions to a BED file
export(cpg_positions, "hg19_all_cpg.bed", format = "BED")

cat("CpG positions saved to 'hg19_all_cpg.bed'\n")
