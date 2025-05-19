#!/bin/bash
#$ -N paleo_checkLibrariesDamage
#$ -S /bin/bash
#$ -l tmem=12G
#$ -l h_vmem=12G
#$ -l h_rt=24:00:00
#$ -wd /SAN/ghlab/pophistory/Alice/paleo_project/logs # one err and out file per sample
#$ -R y # reserve the resources, i.e. stop smaller jobs from getting into the queue while you wait for all the required resources to become available for you

INDIR1="/SAN/ghlab/pophistory/Alice/paleo_project/data/01rawfastq/single/TEMP/01Trimmed_data"
INDIR2="/SAN/ghlab/pophistory/Alice/paleo_project/data/01rawfastq/paired/TEMP/01Trimmed_data"

OUTDIR="/SAN/ghlab/pophistory/Alice/paleo_project/data/02samplesBams"

DamageProfiler="/SAN/ghlab/pophistory/Alice/paleo_project/code/DamageProfiler-1.1-java11.jar"
java="/share/apps/java/bin/java"

## Check if USER/UDR treated with damage profiler --> keep only UDR/USER treated files

for bam in "$INDIR1"/*.filtered.sorted.dedup.bam "$INDIR2"/*.filtered.sorted.dedup.bam; do
    base_name=$(basename "$bam" .filtered.sorted.dedup.bam)
    out_dir="$OUTDIR/damage/output_${base_name}"
    mkdir -p $out_dir
    echo $base_name

    ## 1. Run DamageProfiler (if not done yet):
    if [ ! -f $out_dir/5p_freq_misincorporations.txt ]; then
	$java -jar $DamageProfiler -i "$bam" -o "$out_dir" -sslib
    fi

    ## 2. Extract the C-T info for the 25th first bases:
    # Extract C>T values as a comma-separated string
    ct_values=$(tail -n +5 $out_dir/5p_freq_misincorporations.txt | cut -f2 | paste -sd, -)
    
    # Append to CSV
    echo "$base_name,$ct_values" >> "$OUTDIR/damage/ct_damage_profile.csv"
done

## Remove duplicates
awk 'NR==1{print; next} !seen[$0]++' "$OUTDIR/damage/ct_damage_profile.csv" > TEMP
cat TEMP > "$OUTDIR/damage/ct_damage_profile.csv"
rm TEMP

## If your library has been UDG treated you will expect to see extremely-low to no misincorporations at read ends

Rscript - <<EOF

library(reshape2)
library(ggplot2)

dam <- read.csv("$OUTDIR/damage/ct_damage_profile.csv")
names(dam) <- c("Library", 1:25)
dam <- melt(dam)

pdf(file = "$OUTDIR/damage/plotDamage.pdf", width = 8, height = 4)
ggplot(dam, aes(x = variable, y = value, group = Library, col = Library)) +
  geom_line() + theme_bw() + xlab("position from 5' end") + ylab("C>T misincorporation frequency") +
  theme(legend.position = "none")
dev.off()

EOF

## --> keep only UDR/USER treated files = both values low <0.01

## Untreated ancient DNA shows high C→T at the first position (10–50%) and much lower at the second (1–10%). UDG/USER-treated DNA shows very low C→T (< 1%) at both positions.
## epiPALEOMIX can analyze untreated ancient DNA by leveraging patterns of post-mortem damage and using statistical inference to infer methylation and nucleosome maps, but the results are less precise due to higher background noise from deaminated unmethylated cytosines. USER/UDG treatment is recommended for more accurate methylation inference.

## Additionally, cytosine to uracil deamination that occurs in the ancient DNA molecules due to degradation can be treated using a combination of uracil–DNA–glycosylase (UDG) and endonuclease VIII enzymes, which results in the removal of uracil residues, thus leaving the undamaged parts of the DNA fragments intact (Briggs et al. 2010) as shown in Figure 1.1. Furthermore, partial UDG treatment - by using an E.coli uracil–DNA–glycosylase (UDG-half) - has been developed to retain only deamination on the first and the last two terminal ends of the molecules, thus still enabling the authentication of the recovered molecules (Rohland et al. 2015).

METADATA="/SAN/ghlab/pophistory/Alice/paleo_project/data/00samplesList/metadata_Marchi2022_15samples.tsv"
DAM="$OUTDIR/damage/ct_damage_profile.csv"

## Subset METADATA for all the run_accession for which the second column of DAM < 0.1 (after seeing the damage plot)
OUTPUT="$OUTDIR/damage/filtered_metadata.tsv"

awk -F, 'NR>1 && $2 < 0.1 {print $1}' "$DAM" | 
  awk 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' - "$METADATA" > "$OUTPUT"
