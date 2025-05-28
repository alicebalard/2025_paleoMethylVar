#!/bin/bash
#$ -N paleo_epipaleomix
#$ -S /bin/bash
#$ -l h_vmem=4G,tmem=4G ## NB: tmem value is per core
#$ -pe smp 8 # Request N cores per task 
#$ -t 1-69
#$ -l h_rt=24:00:00
#$ -wd /SAN/ghlab/epigen/Alice/paleo_project/logs # one err and out file per sample
#$ -R y # reserve the resources, i.e. stop smaller jobs from getting into the queue while you wait for all the required resources to become available for you

THREADS=8

# source Python 2.7 to work with epiPALEOMIX
source /share/apps/source_files/python/python-2.7.16.source

DATAPROJ="/SAN/ghlab/epigen/Alice/paleo_project/data"

CODEDIR=$DATAPROJ/../code
BAMDIR="$DATAPROJ/02samplesBams/Antonio2019"
OUTDIR="$DATAPROJ/04methyl/Antonio2019"
## the reference BED file for CpG positions:
BED="$DATAPROJ/03refBed/autosomal_CpGs_hg19_no_dbSNP142.bed"
## Harmonise chromosome names to mach the bam files (done, once)
## sed 's/^ch//' $BED > $BED.harmonised.bed

## https://bitbucket.org/khanghoj/epipaleomix/wiki/Home
## https://www.sciencedirect.com/science/article/pii/S009286742200455X#sec6

## Epipaleomix needs:

## 1. the hg19 reference fasta file: reference/hs37d5.fa.gz
## 2. the bed file: $DATAPROJ/03refBed/autosomal_CpGs_hg19_no_dbSNP142.bed.harmonised.bed
## 3. the library type of the sample (single-stranded/double-stranded): double-stranded for Marchi et al. 2022

## 4. The BAM file:
# Get array of BAM file paths
bam_files=($BAMDIR/*.filtered.sorted.dedup.bam)

# Get array of sample names (basename without extension)
sample_names=()
for bam in "${bam_files[@]}"; do
    sample_names+=( "$(basename "$bam" .filtered.sorted.dedup.bam)" )
done

# SGE_TASK_ID is 1-based; Bash arrays are 0-based, so subtract 1
index=$((SGE_TASK_ID - 1))

# Select the correct BAM path and sample name
MYBAMPATH="${bam_files[$index]}"
SAMPLENAME="${sample_names[$index]}"

# Example: print them
echo "MYBAMPATH: $MYBAMPATH"
echo "SAMPLENAME: $SAMPLENAME"

## index bam file:
samtools index $MYBAMPATH

## All parameters are in a yalm file.
## Prepare yalm file:
# epiPALEOMIX makefile simple > $CODEDIR/makefile.yaml

## Manual modif of yalm file:
YALM0="$CODEDIR/makefile_mytemplate.yaml"

sed "s/SAMPLENAME/$SAMPLENAME/g" $YALM0 | sed "s|MYBAMPATH|$MYBAMPATH|g" - > "$OUTDIR/makefile_$SAMPLENAME.yaml"

## echo "Mean coverage for this sample:" ## needs some G, rather do after on the selected positions
## /share/apps/genomics/bedtools-2.30.0/bin/bedtools coverage -a $BED.harmonised.bed -b $MYBAMPATH | awk '{sum+=$8} END {print sum/NR}'

OUTPUT=$OUTDIR/OUT_makefile_$SAMPLENAME/$SAMPLENAME/${SAMPLENAME}_MethylMap_autosomalCpGshg19noSNP.txt.gz
if [ ! -e "$OUTPUT" ]; then
    echo "Run epiPALEOMIX:" # takes ~2h
    ~/install/epiPALEOMIX/run.py run "$OUTDIR/makefile_$SAMPLENAME.yaml" --max-threads $THREADS --destination $OUTDIR

    echo "Filter output for >=4 reads"
    zcat $OUTPUT | awk 'NR==1 || $4 >= 4' | gzip > $OUTDIR/OUT_makefile_$SAMPLENAME/$SAMPLENAME/${SAMPLENAME}_MethylMap_autosomalCpGshg19noSNP.cov4filtered.txt.gz
else
    echo "epiPALEOMIX already ran on the sample"
fi

## Copy all output in a folder:
mkdir -p $OUTDIR/ALL_MethylMaps
find $OUTDIR -path "$OUTDIR/ALL_MethylMaps" -prune -o \ ## prune to exclude the target dir from the search and -o does the following 
    -type f -name '*_MethylMap_autosomalCpGshg19noSNP.cov4filtered.txt.gz' \
    -exec cp {} $OUTDIR/ALL_MethylMaps/ \;

echo "Analysis finished!"

