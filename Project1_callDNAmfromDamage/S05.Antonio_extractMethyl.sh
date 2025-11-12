#!/bin/bash
#$ -N epipaleomix
#$ -S /bin/bash
#$ -l h_vmem=8G,tmem=8G ## NB: tmem value is per core
#$ -pe smp 6 # Request N cores per task 
#$ -t 1-14
#$ -tc 20
#$ -l h_rt=72:00:00
#$ -wd /SAN/ghlab/epigen/Alice/paleo_project/logs # one err and out file per sample
#$ -R y # reserve the resources, i.e. stop smaller jobs from getting into the queue while you wait for all the required resources to become available for you

THREADS=6
BATCH_SIZE=10 ## Process samples 10 by 10

## Convert SGE_TASK_ID to 0-based
start_index=$(( (SGE_TASK_ID - 1) * BATCH_SIZE ))
end_index=$(( start_index + BATCH_SIZE - 1 ))

# If end_index goes beyond array length, truncate
if [ $end_index -ge ${#bam_files[@]} ]; then
    end_index=$((${#bam_files[@]} - 1))
fi

# source Python 2.7 to work with epiPALEOMIX
source /share/apps/source_files/python/python-2.7.16.source

DATAPROJ="/SAN/ghlab/epigen/Alice/paleo_project/data"

CODEDIR=$DATAPROJ/../code
BAMDIR="$DATAPROJ/02samplesBams/Antonio2019"
OUTDIR="$DATAPROJ/04methyl/Antonio2019/cov1" ## choose coverage here 
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

## 5. All parameters are in a yalm file.
## Prepare yalm file:
# epiPALEOMIX makefile simple > $CODEDIR/makefile.yaml

## echo "Mean coverage for this sample:" ## needs some G, rather do after on the selected positions
## /share/apps/genomics/bedtools-2.30.0/bin/bedtools coverage -a $BED.harmonised.bed -b $MYBAMPATH | awk '{sum+=$8} END {print sum/NR}'

for idx in $(seq $start_index $end_index); do
    MYBAMPATH="${bam_files[$idx]}"
    SAMPLENAME="${sample_names[$idx]}"

    echo "Processing sample: $SAMPLENAME"

    ## The 'OUT_makefileName' folder contains a subfolder for each BAM-file analyzed, including final output files for each analysis conducted in the form of a "BAMName_AnalysisName_BedName.txt.gz" file.
    RAW_OUTPUT=$OUTDIR/OUT_makefile_$SAMPLENAME/$SAMPLENAME/${SAMPLENAME}_MethylMap_autosomalCpGshg19noSNP.txt.gz
    FILTERED_OUTPUT=$OUTDIR/ALL_MethylMaps/${SAMPLENAME}_MethylMap_autosomalCpGshg19noSNP.cov1filtered.txt.gz

    if [ ! -e "$FILTERED_OUTPUT" ]; then
        echo "Make yalm file:"
        YALM0="$CODEDIR/makefile_mytemplate.yaml"
        sed "s/SAMPLENAME/$SAMPLENAME/g" $YALM0 | sed "s|MYBAMPATH|$MYBAMPATH|g" - > "$OUTDIR/makefile_$SAMPLENAME.yaml"

        echo "Index bam:"
        samtools index $MYBAMPATH

        echo "Run epiPALEOMIX:"
        ~/install/epiPALEOMIX/run.py run "$OUTDIR/makefile_$SAMPLENAME.yaml" --max-threads $THREADS --destination $OUTDIR

        echo "Filter output for >=1 reads"
        zcat $RAW_OUTPUT | awk 'NR==1 || $4 >= 1' | gzip > $FILTERED_OUTPUT
    else
        echo "Sample $SAMPLENAME already processed."
    fi

    # Cleanup
    rm -rf $OUTDIR/OUT_makefile_$SAMPLENAME/
    rm -rf $OUTDIR/TEMPORARYFILES_makefile_$SAMPLENAME/
    rm "$OUTDIR/makefile_$SAMPLENAME.yaml"
done

echo "Analysis finished!"

