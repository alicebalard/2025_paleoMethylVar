#!/bin/bash
#$ -N paleo_epipaleomix
#$ -S /bin/bash
#$ -l h_vmem=4G,tmem=4G ## NB: tmem value is per core
#$ -pe smp 8 # Request N cores per task 
#$ -t 1-15
#$ -l h_rt=100:00:00
#$ -wd /SAN/ghlab/pophistory/Alice/paleo_project/logs # one err and out file per sample
#$ -R y # reserve the resources, i.e. stop smaller jobs from getting into the queue while you wait for all the required resources to become available for you

THREADS=8

# source Python 2.7 to work with epiPALEOMIX
source /share/apps/source_files/python/python-2.7.16.source

DATAPROJ="/SAN/ghlab/pophistory/Alice/paleo_project/data"

CODEDIR=$DATAPROJ/../code
BAMDIR="$DATAPROJ/02samplesBams"
OUTDIR="$DATAPROJ/04methyl"
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
# Create the files to loop over if it does not exist:
if [ ! -e "$BAMDIR/list_of_bams.tmp" ]; then
    ls -1 $BAMDIR/*.sorted.dedup.bam > $BAMDIR/list_of_bams.tmp
    xargs -n1 basename < $BAMDIR/list_of_bams.tmp | cut -d. -f1 > $BAMDIR/list_of_sample_names.tmp
fi

## Select the correct line of list of files at each iteration
MYBAMPATH=$(sed -n "${SGE_TASK_ID}p" $BAMDIR/list_of_bams.tmp)
SAMPLENAME=$(sed -n "${SGE_TASK_ID}p" $BAMDIR/list_of_sample_names.tmp)

## All parameters are in a yalm file.
## Prepare yalm file:
# epiPALEOMIX makefile simple > $CODEDIR/makefile.yaml

## Manual modif of yalm file:
YALM0="$CODEDIR/makefile_template_Marchi.yaml"

sed "s/SAMPLENAME/$SAMPLENAME/g" $YALM0 | sed "s|MYBAMPATH|$MYBAMPATH|g" - > "$OUTDIR/makefile_Marchi_$SAMPLENAME.yaml"

echo "Mean coverage for this sample:" ## needs some G
/share/apps/genomics/bedtools-2.30.0/bin/bedtools coverage -a $BED.harmonised.bed -b $MYBAMPATH | awk '{sum+=$8} END {print sum/NR}'

echo "Run epiPALEOMIX:" # takes ~2h
~/install/epiPALEOMIX/run.py run "$OUTDIR/makefile_Marchi_$SAMPLENAME.yaml" --max-threads $THREADS --destination $OUTDIR

OUTPUT=$OUTDIR/OUT_makefile_Marchi_$SAMPLENAME/$SAMPLE_NAME/${SAMPLENAME}_MethylMap_autosomalCpGshg19noSNP.txt.gz

echo "Filter output for >=4 reads"

zcat $OUTPUT | awk 'NR==1 || $4 >= 4' | gzip > $OUTDIR/OUT_makefile_Marchi_$SAMPLENAME/$SAMPLE_NAME/${SAMPLENAME}_MethylMap_autosomalCpGshg19noSNP.cov4filtered.txt.gz

echo "Analysis finished!"



