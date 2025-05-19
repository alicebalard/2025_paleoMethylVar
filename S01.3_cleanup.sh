#!/bin/bash

# Usage: ./S01.3_cleanup.sh /path/to/your/working/directory

# Check if a directory argument is provided
if [[ -z "$1" ]]; then
    echo "Usage: $0 /path/to/working/directory"
    exit 1
fi

DIR="$1"

# Ensure the provided directory exists
if [[ ! -d "$DIR" ]]; then
    echo "Error: Directory '$DIR' does not exist."
    exit 1
fi

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
