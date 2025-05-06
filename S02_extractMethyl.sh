# Step 4: CpG Position Handling (R Code)
echo "Generating CpG positions..."
Rscript - <<EOF
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


###################################
## Bed file to use for epiPaleomix:

## "We called all CG dinucleotide autosomal positions (n = 26,752,702) from the human (hg19) reference genome using the R Bioconductor package “BSgenome.Hsapiens.UCSC.hg19” (Pagès, 2019) and stored these in a BED file. We then filtered these by removing any positions overlapping with SNP positions from dbSNP 142 (Sherry et al., 2001). Our aim here was to avoid confounding between methylation signals and real variants at CpG positions. There remained 13,270,411 autosomal CpG positions in the reference genome."

# Download dbSNP142 (requires appropriate download URL)
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/BED/bed_chr_1.bed.gz
# Repeat for all chromosomes and concatenate

# Use bedtools to filter
bedtools intersect -a hg19_autosomal_cpg.bed \
                   -b dbSNP142_combined.bed \
                   -v > filtered_cpg_no_snp.bed

