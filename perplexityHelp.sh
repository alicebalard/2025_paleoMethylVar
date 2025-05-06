#!/bin/bash

# Configuration
REFERENCE_HS37D5="hs37d5.fa"
REFERENCE_HG19="hg19.fa"
THREADS=8
ENA_ACCESSION_LIST="Table_S1_accessions.txt"
OUTPUT_DIR="analysis_results"

# Create output directory
mkdir -p $OUTPUT_DIR

# Step 1: Data Download from ENA
echo "Downloading data from ENA..."
while read -r acc; do
    # Download FASTQ/BAM files
    if [[ $acc =~ ^ERR ]]; then
        fastq-dump --split-files --gzip $acc
    elif [[ $acc =~ ^DRR ]]; then
        ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
            era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/${acc:0:6}/$acc/ $OUTPUT_DIR/
    fi
done < $ENA_ACCESSION_LIST

# Step 2: Remapping with bwa aln
process_sample() {
    local input=$1
    local sample_name=$(basename $input | cut -d. -f1)
    
    echo "Processing $sample_name..."
    
    # Convert BAM to FASTQ if needed
    if [[ $input == *.bam ]]; then
        samtools fastq $input > $OUTPUT_DIR/${sample_name}.fastq
        input=$OUTPUT_DIR/${sample_name}.fastq
    fi

    # Alignment
    bwa aln -t $THREADS -l 16500 -n 0.01 -o 2 $REFERENCE_HS37D5 $input > $OUTPUT_DIR/${sample_name}.sai
    bwa samse $REFERENCE_HS37D5 $OUTPUT_DIR/${sample_name}.sai $input > $OUTPUT_DIR/${sample_name}.sam

    # Convert to BAM
    samtools view -@ $THREADS -Sb $OUTPUT_DIR/${sample_name}.sam > $OUTPUT_DIR/${sample_name}.bam

    # Filtering
    samtools view -h $OUTPUT_DIR/${sample_name}.bam | \
        awk 'length($10) >= 35 || $1 ~ /^@/' | \
        samtools view -b -q 30 - | \
        samtools view -h - | \
        awk 'BEGIN {mismatch_field=12} 
             $1 ~ /^@/ || ($mismatch_field !~ /NM:i:/) {print; next} 
             {split($mismatch_field, a, ":"); mismatches=a[3]; len=length($10); 
              if (mismatches/len <= 0.1) print}' | \
        samtools view -Sb - > $OUTPUT_DIR/${sample_name}_filtered.bam

    # Sort and index
    samtools sort -@ $THREADS $OUTPUT_DIR/${sample_name}_filtered.bam -o $OUTPUT_DIR/${sample_name}_final.bam
    samtools index $OUTPUT_DIR/${sample_name}_final.bam
}

# Process all samples
for file in *.fastq.gz *.bam; do
    process_sample "$file"
done

# Step 3: PMD Analysis
echo "Running PMDtools..."
for bam in $OUTPUT_DIR/*_final.bam; do
    pmdtools --deamination > $OUTPUT_DIR/$(basename $bam .bam)_pmd.txt
done

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

# Step 5: Methylation Analysis with epiPALEOMIX
echo "Running epiPALEOMIX..."
for bam in $OUTPUT_DIR/*_final.bam; do
    sample=$(basename $bam .bam)
    epiPALEOMIX analyze \
        --bam $bam \
        --reference $REFERENCE_HG19 \
        --cpg-bed hg19_cpg_filtered.bed \
        --library-type double-stranded \
        --min-coverage 2 \
        --output $OUTPUT_DIR/${sample}_methylation.txt
done

# Step 6: Generate Combined Methylation Matrix (R Code)
echo "Generating methylation matrix..."
Rscript - <<EOF
library(dplyr)
library(tidyr)

file_list <- list.files(pattern="*_methylation.txt", full.names=TRUE)
all_data <- lapply(file_list, function(f) {
    read.table(f, header=TRUE) %>%
    filter(coverage >= 4) %>%
    select(chr, pos, methylation_score) %>%
    rename(!!basename(f) := methylation_score)
})

combined <- all_data %>% reduce(full_join, by=c("chr", "pos"))
write.csv(combined, "methylation_matrix.csv", na="NA", row.names=FALSE)

# Calculate average methylation scores
average_scores <- combined %>%
    gather(sample, score, -chr, -pos) %>%
    group_by(chr, pos) %>%
    summarise(mean_score = mean(score, na.rm=TRUE))

write.csv(average_scores, "average_methylation_scores.csv", row.names=FALSE)
EOF

# Step 7: Region-based Analysis (R Code)
echo "Analyzing CGI regions..."
Rscript - <<EOF
library(GenomicRanges)
library(rtracklayer)

# Import regions
regions <- list(
    CGI = import("cpgIslandExt.hg19.bed"),
    Shores5 = import("cgi_shores5.bed"),
    Shores3 = import("cgi_shores3.bed"),
    Shelves5 = import("cgi_shelves5.bed"),
    Shelves3 = import("cgi_shelves3.bed")
)

# Calculate region averages
methylation <- read.csv("average_methylation_scores.csv")
methyl_gr <- makeGRangesFromDataFrame(methylation, 
    seqnames.field="chr", start.field="pos", end.field="pos")

results <- lapply(regions, function(region) {
    overlaps <- findOverlaps(methyl_gr, region)
    methylation[queryHits(overlaps),] %>%
        group_by(seqnames) %>%
        summarise(mean_methylation = mean(mean_score, na.rm=TRUE))
})

# Save results
mapply(function(data, name) {
    write.csv(data, paste0("methylation_", name, ".csv"), row.names=FALSE)
}, results, names(regions))
EOF

Key Components Explained
Data Download

Handles both FASTQ (via fastq-dump) and BAM (via Aspera) downloads

Uses accession list from Table S1

Alignment & Filtering

Implements strict filtering (≥35bp, MAPQ≥30, ≤10% mismatches)

Maintains proper BAM sorting/indexing

PMD Analysis

Uses pmdtools for deamination pattern verification

Generates PMD profiles for damage assessment

CpG Processing

Uses BSgenome.Hsapiens.UCSC.hg19 to identify CpGs

Filters against dbSNP142 variants

Annotates CGI regions (islands, shores, shelves)

Methylation Analysis

Runs epiPALEOMIX with coverage filters (≥2×, ≥4 reads per CpG)

Generates sample-specific methylation scores

Data Integration

Creates merged methylation matrix with NA handling

Calculates position-specific averages

Region-based Analysis

Computes average methylation for different CGI regions

Generates separate reports for islands/shores/shelves
