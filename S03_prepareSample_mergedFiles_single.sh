#!/bin/bash
#$ -N paleo_prepSample_single_check_damage
#$ -S /bin/bash
#$ -l tmem=2G
#$ -l h_vmem=2G
#$ -l h_rt=240:00:00
#$ -wd /SAN/ghlab/pophistory/Alice/paleo_project/logs # one err and out file per sample
#$ -R y # reserve the resources, i.e. stop smaller jobs from getting into the queue while you wait for all the required resources to become available for you

INDIR="/SAN/ghlab/pophistory/Alice/paleo_project/data/01rawfastq/single/TEMP/01Trimmed_data"
OUTDIR="/SAN/ghlab/pophistory/Alice/paleo_project/data/02samplesBams"
DamageProfiler="/SAN/ghlab/pophistory/Alice/paleo_project/code/DamageProfiler-1.1-java11.jar"
java="/share/apps/java/bin/java"

## Check if USER/UDR treated with damage profiler --> keep only UDR/USER treated files
output_csv="$INDIR/damage_report.csv"
echo "BAM_Name,Pos1_CT,Pos2_CT" > "$output_csv"

for bam in $INDIR/*.fq.gz.sorted.dedup.bam; do
    base_name=$(basename "$bam" .fq.gz.sorted.dedup.bam)
    out_dir="$INDIR/output_${base_name}"
    
    # Run DamageProfiler
    $java -jar $DamageProfiler -i "$bam" -o "$out_dir" -sslib
    
#    cd $out_dir 
#    # Extract values (first few rows are just text)
#    pos1=$(awk 'NR==5{print $2}' "$out_dir/5p_freq_misincorporations.txt")
#    pos2=$(awk 'NR==6{print $2}' "$out_dir/5p_freq_misincorporations.txt")
    
    # Append to CSV
#    echo "${base_name},${pos1},${pos2}" >> "$output_csv"
done

## If your library has been UDG treated you will expect to see extremely-low to no misincorporations at read ends

## --> keep only UDR/USER treated files = both values low <0.01

#awk -F, 'NR==1 || ($2 > 0.05 && $3 < 0.02)' damage_summary.csv > partial_udg_samples.csv
#Enough terminal damage is present (Pos1 > 5%)
#Internal damage has been removed (Pos2 is low)

## Untreated ancient DNA shows high C→T at the first position (10–50%) and much lower at the second (1–10%). UDG/USER-treated DNA shows very low C→T (< 1%) at both positions.
## epiPALEOMIX can analyze untreated ancient DNA by leveraging patterns of post-mortem damage and using statistical inference to infer methylation and nucleosome maps, but the results are less precise due to higher background noise from deaminated unmethylated cytosines. USER/UDG treatment is recommended for more accurate methylation inference.

## Merge the bams by sample
/SAN/ghlab/pophistory/Alice/paleo_project/data/01rawfastq/single/TEMP/01Trimmed_data/ERR9237189_trimmed.fq.gz.sorted.bam

#!/bin/bash

# Paths (adjust as needed)
BAM_DIR="/SAN/ghlab/pophistory/Alice/paleo_project/data/01rawfastq/single/TEMP/01Trimmed_data"
METADATA="metadata.tsv"  # Your metadata file path
OUTPUT_DIR="./merged_bams"
mkdir -p "$OUTPUT_DIR"

# Skip header, read run_accession and sample_accession, group runs by sample
declare -A sample_runs

tail -n +2 "$METADATA" | while IFS=$'\t' read -r run_accession experiment_accession sample_alias study_accession fastq_ftp sample_accession; do
    sample_runs["$sample_accession"]+="$run_accession "
done

# For each sample, merge BAMs
for sample in "${!sample_runs[@]}"; do
    echo "Processing sample: $sample"
    runs=${sample_runs[$sample]}
    
    bam_files=()
    for run in $runs; do
        bam_file="$BAM_DIR/${run}_trimmed.fq.gz.sorted.bam"
        if [[ -f "$bam_file" ]]; then
            bam_files+=("$bam_file")
        else
            echo "Warning: BAM file not found: $bam_file"
        fi
    done

    if [ ${#bam_files[@]} -eq 0 ]; then
        echo "No BAM files found for sample $sample, skipping."
        continue
    fi

    # Merge BAM files using samtools
    merged_bam="$OUTPUT_DIR/${sample}.merged.bam"
    echo "Merging ${#bam_files[@]} BAM files into $merged_bam"
    samtools merge -@ 4 "$merged_bam" "${bam_files[@]}"
    
    # Optional: index merged BAM
    samtools index "$merged_bam"
done

## One last duplicate removal

