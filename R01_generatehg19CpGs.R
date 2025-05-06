# We called all CG dinucleotide autosomal positions (n = 26,752,702) from the human (hg19) reference genome using the R Bioconductor package “BSgenome.Hsapiens.UCSC.hg19” (Pagès, 2019) and stored these in a BED file. We then filtered these by removing any positions overlapping with SNP positions from dbSNP 142 (Sherry et al., 2001). Our aim here was to avoid confounding between methylation signals and real variants at CpG positions. There remained 13,270,411 autosomal CpG positions in the reference genome.
# 
# We downloaded CpG island (CGI) positions for hg19 from the UCSC Genome Browser (Karolchik et al., 2004). We termed 2 kb sequences flanking CpG islands “shores” (upstream regions “shores5” and downstream regions “shores3”), 2 kb sequences flanking the shores “shelves” (upstream regions “shelves5” and downstream regions “shelves3”), and distal sites outside the CpG island regions as “open sea,” following (Hanghøj et al., 2016).



## Step 4: CpG Position Handling (R Code)
library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)

# Get all CpG positions
cpg_positions <- which(vmatchPattern("CG", BSgenome.Hsapiens.UCSC.hg19\$autosomes) != -1)
cpg_bed <- GRanges(names(cpg_positions), IRanges(start=cpg_posions, width=2))
export.bed(cpg_bed, "hg19_cpg_all.bed")

# Filter against dbSNP
dbsnp <- import("dbsnp142.hg19.bed")
filtered_cpg <- subsetByOverlaps(cpg_bed, dbsnp, invert=TRUE)
export.bed(filtered_cpg, "hg19_cpg_filtered.bed")

# Annotate regions
cgi <- import("cpgIslandExt.hg19.bed")
flank <- 2000
shores5 <- flank(cgi, flank, start=TRUE)
shores3 <- flank(cgi, flank, start=FALSE)
shelves5 <- flank(shores5, flank, start=TRUE)
shelves3 <- flank(shores3, flank, start=FALSE)

export.bed(shores5, "cgi_shores5.bed")
export.bed(shores3, "cgi_shores3.bed")
export.bed(shelves5, "cgi_shelves5.bed")
export.bed(shelves3, "cgi_shelves3.bed")
EOF



########### previous
## We called all CG dinucleotide autosomal positions (n = 26,752,702) from the human (hg19) reference genome using the R Bioconductor package “BSgenome.Hsapiens.UCSC.hg19” (Pagès, 2019) and stored these in a BED file. 

## Install required Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg19", "GenomicRanges", "rtracklayer"))

library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(rtracklayer)

# Get autosomal chromosomes
autosomal_chr <- paste0("chr", 1:22)

# Find CG positions chromosome-by-chromosome
cgranges <- GRanges()
for (chr in autosomal_chr) {
    cat("Processing", chr, "\n")
    chr_seq <- BSgenome.Hsapiens.UCSC.hg19[[chr]]
    cg_matches <- matchPattern("CG", chr_seq)
    cgranges <- c(cgranges, GRanges(
        seqnames = chr,
        ranges = IRanges(start = start(cg_matches), width = 2),
        strand = "*"
    ))
}

# Remove duplicates and sort
cgranges <- unique(sort(cgranges))

# Export to BED file
export(cgranges, "hg19_autosomal_cpg.bed", format = "BED")
