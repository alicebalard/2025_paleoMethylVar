#!/usr/bin/env bash
set -euo pipefail

# === 0. Go to the folder with your trimmed FASTQs ===
# Edit this to your actual path (WSL-style) or comment it out if you run the script from there already.
# cd /mnt/c/Users/hamda/Desktop/UGI_thingy/3rdyear/genome/quality

# === 1. Make an adapter file with all over-represented sequences ===
cat > overrep_adapters.fa << 'ADP'
>overrep_CGAGGAT_nohit
CGAGGATTCGAAAAGGTGAACCGAGCCAGAC

>overrep_GGCTTCTT_nohit
GGCTTCTTCTTGCTCTTAGGTAGTAGGTTC

>overrep_10658_1
TTGAGATTTTTGTTAAAGGAGATGTTTGTGGTGTGTATATGTTTTATGTT

>overrep_10658_2
ATGTTTGTGGTGTGTATATGTTTTATGTTGGTTGGTTGAGTATTTTGGTG

>overrep_TruSeq_GCCTAC
AAGAGCACACGTCTGAACTCCAGTCACGCCTACGATCTCGTATGCCGTCT

>overrep_TruSeq_TAGGCC
AAGAGCACACGTCTGAACTCCAGTCACTAGGCCGATCTCGTATGCCGTCT

>overrep_TruSeq_GACGATTAT_1
AAGAGCACACGTCTGAACTCCAGTCACGACGATTATCTCGTATGCCGTCT

>overrep_TruSeq_GACGATTAT_2
AGAGCACACGTCTGAACTCCAGTCACGACGATTATCTCGTATGCCGTCTT

>overrep_TruSeq_GACGATTAT_3
AAGAGCACACGTCTGAACTCCAGTCACGACGATTATCGCGTATGCCGTCT

>overrep_TruSeq_GACGATTAT_4
GAAGAGCACACGTCTGAACTCCAGTCACGACGATTATCTCGTATGCCGTC

>overrep_TruSeq_GCAGAAG
AAGAGCACACGTCTGAACTCCAGTCACGCAGAAGATCTCGTATGCCGTCT
ADP

echo "Written overrep_adapters.fa with common over-represented sequences"

# === 2. Loop over all trimmed FastQ files and trim these sequences ===
shopt -s nullglob

for fq in *_bisilfite.trim.fastq.gz; do
  [ -e "$fq" ] || continue

  base=${fq%.fastq.gz}          # e.g. Vac_210_Molar_bisilfite.trim
  out="${base}.overrep_trim.fastq.gz"

  echo "Trimming over-represented sequences from: $fq"
  cutadapt \
    -j 8 \
    -b file:overrep_adapters.fa \
    -m 25 \
    -o "$out" \
    "$fq"

  echo "  -> wrote $out"
done

echo "Done trimming over-represented sequences."
