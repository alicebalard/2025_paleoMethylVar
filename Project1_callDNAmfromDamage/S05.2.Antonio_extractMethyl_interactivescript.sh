#!/bin/bash

THREADS=8

# source Python 2.7 for epiPALEOMIX
source /share/apps/source_files/python/python-2.7.16.source

DATAPROJ="/SAN/ghlab/epigen/Alice/paleo_project/data"

CODEDIR=$DATAPROJ/../code
BAMDIR="$DATAPROJ/02samplesBams/Antonio2019"
OUTDIR="$DATAPROJ/04methyl/Antonio2019/noCovFilter"
BED="$DATAPROJ/03refBed/autosomal_CpGs_hg19_no_dbSNP142.bed"

mkdir -p $OUTDIR/ALL_MethylMaps

# Get array of BAM files and sample names
bam_files=($BAMDIR/*.filtered.sorted.dedup.bam)
sample_names=()
for bam in "${bam_files[@]}"; do
    sample_names+=( "$(basename "$bam" .filtered.sorted.dedup.bam)" )
done

echo "Starting epiPALEOMIX loop for ${#sample_names[@]} samples."
echo "Job started at: $(date)"
job_start=$(date +%s)

# Loop over all samples
for i in "${!sample_names[@]}"; do
    sample=${sample_names[$i]}
    bam=${bam_files[$i]}
    output="$OUTDIR/ALL_MethylMaps/${sample}_MethylMap_autosomalCpGshg19noSNP.cov1filtered.txt.gz"

    echo "[$(date)] Starting sample: $sample"
    sample_start=$(date +%s)

    if [ ! -e "$output" ]; then
        echo "Make YAML file..."
        YALM0="$CODEDIR/makefile_mytemplate.yaml"
        sed "s/SAMPLENAME/$sample/g" $YALM0 | sed "s|MYBAMPATH|$bam|g" - > "$OUTDIR/makefile_$sample.yaml"

        echo "Index BAM..."
        samtools index $bam

        echo "Run epiPALEOMIX..."
        ~/install/epiPALEOMIX/run.py run "$OUTDIR/makefile_$sample.yaml" \
            --max-threads $THREADS \
            --destination $OUTDIR

        echo "Filter output for >=1 reads..."
        zcat "$output" | awk 'NR==1 || $4 >= 1' | gzip > \
            "$OUTDIR/OUT_makefile_$sample/$sample/${sample}_MethylMap_autosomalCpGshg19noSNP.cov1filtered.txt.gz"
    else
        echo "epiPALEOMIX already ran on $sample"
    fi

    echo "Collect results..."
    find "$OUTDIR" -path "$OUTDIR/ALL_MethylMaps" -prune -o \
         -type f -name '*_MethylMap_autosomalCpGshg19noSNP.cov1filtered.txt.gz' \
         -exec cp -uv {} "$OUTDIR/ALL_MethylMaps/" \;

    echo "Clean up temporary files..."
    rm -rf "$OUTDIR/OUT_makefile_$sample/"
    rm -rf "$OUTDIR/TEMPORARYFILES_makefile_$sample/"
    rm -f  "$OUTDIR/makefile_$sample.yaml"

    sample_end=$(date +%s)
    runtime=$((sample_end - sample_start))
    echo "[$(date)] Finished sample: $sample (took $((runtime/60)) min $((runtime%60)) sec)"
done

job_end=$(date +%s)
total_runtime=$((job_end - job_start))
echo "All samples finished at: $(date)"
echo "Total runtime: $((total_runtime/3600)) h $(((total_runtime%3600)/60)) min $((total_runtime%60)) sec"
