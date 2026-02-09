mkdir -p ~/refs/hg38
cd ~/refs/hg38

# Download reference FASTA (UCSC version)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Prepare bisulfite genome
bismark_genome_preparation --verbose ~/refs/hg38b