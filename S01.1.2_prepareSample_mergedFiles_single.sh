#!/bin/bash
#$ -N paleo_prepSample_single_merge
#$ -S /bin/bash
#$ -l tmem=8G
#$ -l h_vmem=8G
#$ -pe smp 8 # Request N cores per task 
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
    
    cd $out_dir 
    # Extract values (first few rows are just text)
    pos1=$(awk 'NR==5{print $2}' "$out_dir/5p_freq_misincorporations.txt")
    pos2=$(awk 'NR==6{print $2}' "$out_dir/5p_freq_misincorporations.txt")
    
    # Append to CSV
    echo "${base_name},${pos1},${pos2}" >> "$output_csv"
done

## If your library has been UDG treated you will expect to see extremely-low to no misincorporations at read ends

## --> keep only UDR/USER treated files = bos values <0.01

## Merge the bams by sample

## One last duplicate removal

