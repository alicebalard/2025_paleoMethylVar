
#################################################################################
## Dl data reviewed in https://onlinelibrary.wiley.com/doi/full/10.1111/eva.13743 

### 15 ancient genomes from Marchi et al. 2022
### https://www.sciencedirect.com/science/article/pii/S009286742200455X#app2
### ENA PRJEB50857

cd /SAN/ghlab/pophistory/Alice/paleo_project/data/00samplesList

wget -O metadata_Marchi2022_15samples.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB50857&result=read_run&fields=experiment_accession,run_accession,sample_alias,study_accession,fastq_ftp,sample_accession"

# Extract all FTP links (including paired-end files)
## tail -n +2 filename prints the file from line 2 to the end.
tail -n +2 metadata_Marchi2022_15samples.tsv | cut -f5 | tr ';' '\n' | sed 's/^/ftp:\/\//' > ftp_links_Marchi2022_15samples.tsv

# Download all files sequentially
cd /SAN/ghlab/pophistory/Alice/paleo_project/data/01rawfastq

while read url; do
    wget --continue --no-clobber "$url"
done < /SAN/ghlab/pophistory/Alice/paleo_project/data/00samplesList/ftp_links_Marchi2022_15samples.tsv

## Check the number of files downloaded:
## cat /SAN/ghlab/pophistory/Alice/paleo_project/data/00samplesList/ftp_links_Marchi2022_15samples.tsv | wc -l ## 878
## ll /SAN/ghlab/pophistory/Alice/paleo_project/data/01rawfastq/*fastq.gz | wc -l # 878

# Organise file types per library
mkdir -p single paired

# Classify based on number of files per sample
awk -F';' '{print NF, $0}' metadata_Marchi2022_15samples.tsv | while read count line; do ## add 1 or 2 at line start, 2 if one ";" cut the line in 2
    if [ $count -eq 1 ]; then
	mv $(echo "$line" | cut -f5 | cut -d';' -f1 | sed 's|.*/||') single/
    elif [ $count -eq 2 ]; then
	mv $(echo "$line" | cut -f5 | cut -d';' -f1 | sed 's|.*/||') paired/
	mv $(echo "$line" | cut -f5 | cut -d';' -f2 | sed 's|.*/||') paired/
    fi
done

# Organise files per sample





