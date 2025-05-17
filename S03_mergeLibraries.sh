#!/bin/bash
#$ -N paleo_mergeLibs
#$ -S /bin/bash
#$ -l tmem=100G
#$ -l h_vmem=100G
#$ -l h_rt=48:00:00
#$ -wd /SAN/ghlab/pophistory/Alice/paleo_project/logs # one err and out file per sample
#$ -R y # reserve the resources, i.e. stop smaller jobs from getting into the queue while you wait for all the required resources to become available for you

## Libraries are in:
#/SAN/ghlab/pophistory/Alice/paleo_project/data/01rawfastq/single/TEMP/01Trimmed_data/*.filtered.sorted.dedup.bam
#/SAN/ghlab/pophistory/Alice/paleo_project/data/01rawfastq/paired/TEMP/01Trimmed_data/*.filtered.sorted.dedup.bam

## 1. Merge the bams libraries by sample
## 2. Remove those non UDG treated (see previous script for filtered metadata)
## 3. Merge by sample_alias

# Set paths
METADATA="/SAN/ghlab/pophistory/Alice/paleo_project/data/02samplesBams/damage/filtered_metadata.tsv"
OUTDIR="/SAN/ghlab/pophistory/Alice/paleo_project/data/02samplesBams/"
SINGLE_DIR="/SAN/ghlab/pophistory/Alice/paleo_project/data/01rawfastq/single/TEMP/01Trimmed_data"
PAIRED_DIR="/SAN/ghlab/pophistory/Alice/paleo_project/data/01rawfastq/paired/TEMP/01Trimmed_data"

mkdir -p "$OUTDIR"

# Extract unique sample_aliases from metadata (skip header)
awk 'NR>1 {print $3}' "$METADATA" | sort | uniq | while read SAMPLE; do
    # For each sample_alias, find all run_accession for that sample
    RUNS=$(awk -v s="$SAMPLE" 'NR>1 && $3==s {print $1}' "$METADATA")
    # For each run_accession, find the corresponding BAM files
    BAMFILES=()
    for RUN in $RUNS; do
        for BAM in "$SINGLE_DIR/${RUN}.filtered.sorted.dedup.bam" "$PAIRED_DIR/${RUN}.filtered.sorted.dedup.bam"; do
            [ -f "$BAM" ] && BAMFILES+=("$BAM")
        done
    done
    # Only merge if we have more than zero BAM files
    if [ ${#BAMFILES[@]} -gt 0 ]; then
        echo "Merging ${#BAMFILES[@]} BAMs for sample $SAMPLE..."
        samtools merge -f "$OUTDIR/${SAMPLE}.bam" "${BAMFILES[@]}"

	## Deduplicate the merged bams
	echo "Sort and index"
	samtools sort -o "$OUTDIR/${SAMPLE}.sorted.bam" "$OUTDIR/${SAMPLE}.bam"
	samtools index "$OUTDIR/${SAMPLE}.sorted.bam"

	echo "Duplicates removal:"
	java="/share/apps/java/bin/java"
	"$java" -jar -Xmx3G -Xms3G /share/apps/genomics/picard-2.20.3/bin/picard.jar MarkDuplicates \
		I="$OUTDIR/${SAMPLE}.sorted.bam" \
		O="$OUTDIR/${SAMPLE}.sorted.dedup.bam" \
		VALIDATION_STRINGENCY=LENIENT \
		REMOVE_DUPLICATES=TRUE \
		M="$OUTDIR/${SAMPLE}.sorted.bam_dup_metrics.txt"

	echo "Index bam file:"
	samtools index "$OUTDIR/${SAMPLE}.sorted.dedup.bam"

	echo "Remove intermediate files:"
	rm "$OUTDIR/${SAMPLE}.bam"
	rm "$OUTDIR/${SAMPLE}.sorted.bam"
	rm "$OUTDIR/${SAMPLE}.sorted.bam.bai"

    else
        echo "No BAMs found for sample $SAMPLE, skipping."
    fi
done
