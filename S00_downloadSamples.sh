####################
### Data downloaded:
### 15 ancient genomes from Marchi et al. 2022
### https://www.sciencedirect.com/science/article/pii/S009286742200455X#app2
### ENA PRJEB50857

### 123 ancient + modern genomes from Antonio et al. 2019 (partial-UDR)
### https://www.science.org/doi/10.1126/science.aay6826#bibliography
### PRJEB32566ENA

DATADIR="/SAN/ghlab/epigen/Alice/paleo_project/data/"
SAMPLELIST_DIR="$DATADIR/00samplesList"
mkdir -p $SAMPLELIST_DIR
cd $SAMPLELIST_DIR

## 1. DL metadata:
wget -O $SAMPLELIST_DIR/metadata_Marchi2022_15samples.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB50857&result=read_run&fields=experiment_accession,run_accession,sample_alias,study_accession,fastq_ftp,sample_accession"

wget -O $SAMPLELIST_DIR/metadata_Antonio2019_127samples.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB32566&result=read_run&fields=experiment_accession,run_accession,sample_alias,study_accession,fastq_ftp,sample_accession"

## 2. Extract all FTP links (including paired-end files)
## tail -n +2 filename prints the file from line 2 to the end.
tail -n +2 $SAMPLELIST_DIR/metadata_Marchi2022_15samples.tsv | cut -f5 | tr ';' '\n' | sed 's/^/ftp:\/\//' > $SAMPLELIST_DIR/ftp_links_Marchi2022_15samples.tsv

tail -n +2 $SAMPLELIST_DIR/metadata_Antonio2019_127samples.tsv | cut -f5 | tr ';' '\n' | sed 's/^/ftp:\/\//' > $SAMPLELIST_DIR/ftp_links_Antonio2019_127samples.tsv

# Download all files sequentially
SAMPLEFASTQ_DIR="$DATADIR/01rawfastq/Marchi2022"
cd $SAMPLEFASTQ_DIR
while read url; do
    wget --continue --no-clobber "$url"
done < $SAMPLELIST_DIR/ftp_links_Marchi2022_15samples.tsv

SAMPLEFASTQ_DIR="$DATADIR/01rawfastq/Antonio2019"
mkdir -p $SAMPLEFASTQ_DIR
cd $SAMPLEFASTQ_DIR
while read url; do
    wget --continue --no-clobber "$url"
done < $SAMPLELIST_DIR/ftp_links_Antonio2019_127samples.tsv

## Check the number of files downloaded:
## cat $SAMPLELIST_DIR/ftp_links_Marchi2022_15samples.tsv | wc -l ## 878
## ll $SAMPLEFASTQ_DIR/*fastq.gz | wc -l # 878

# Organise file types per library
mkdir -p $SAMPLEFASTQ_DIR/single $SAMPLEFASTQ_DIR/paired

# Classify based on number of files per sample
cd $SAMPLEFASTQ_DIR
awk -F';' '{print NF, $0}' $SAMPLELIST_DIR/metadata_Marchi2022_15samples.tsv | while read count line; do ## add 1 or 2 at line start, 2 if one ";" cut the line in 2
    if [ $count -eq 1 ]; then
	mv $(echo "$line" | cut -f5 | cut -d';' -f1 | sed 's|.*/||') $SAMPLEFASTQ_DIR/single/
    elif [ $count -eq 2 ]; then
	mv $(echo "$line" | cut -f5 | cut -d';' -f1 | sed 's|.*/||') $SAMPLEFASTQ_DIR/paired/
	mv $(echo "$line" | cut -f5 | cut -d';' -f2 | sed 's|.*/||') $SAMPLEFASTQ_DIR/paired/
    fi
done

## For Antonio 2019, only single end
