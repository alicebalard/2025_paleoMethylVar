#!/bin/bash

# Set the working directory
DIR="/SAN/ghlab/pophistory/Alice/paleo_project/data/01rawfastq/single/TEMP/01Trimmed_data"

cd "$DIR" || exit 1

# Loop over all .filtered.sorted.dedup.bam files (get SAMPLENAME)
for bam in *.filtered.sorted.dedup.bam; do
    # Skip if no such files
    [[ -e "$bam" ]] || continue

    # Extract SAMPLENAME (everything before .filtered.sorted.dedup.bam)
    samplename="${bam%.filtered.sorted.dedup.bam}"

    echo "Keeping: $bam"
    echo "Removing all: ${samplename}_* and ${samplename}.* except $bam"

    # Find all files matching SAMPLENAME_* or SAMPLENAME.*
    for f in ${samplename}_* ${samplename}.*; do
        # Skip the file to keep
        [[ "$f" == "$bam" ]] && continue
        # Skip if file does not exist (globbing may return pattern if no match)
        [[ -e "$f" ]] || continue
        echo "Deleting $f"
        rm -f -- "$f"
    done
done


