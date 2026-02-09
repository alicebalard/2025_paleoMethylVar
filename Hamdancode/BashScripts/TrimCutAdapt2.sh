#!/usr/bin/env bash
set -euo pipefail

# ==== CONFIG ====
THREADS=8
MIN_PCT=0.20      # only trim overrepresented sequences with >= this % in FastQC table
MIN_LEN=30        # throw away reads shorter than this after trimming
DATA_DIR="$(pwd)" # run the script from the folder with your FASTQs

# Optional: also trim homopolymer tails common on some runs
TRIM_HOMOPOLYMERS=true
HOMOPOLYMER_LEN=18   # e.g. G{18} and A{18}

# ==== PREP ====
mkdir -p qc adapters trimmed logs

# Utility: extract Overrepresented sequences from a FastQC 'fastqc_data.txt'
extract_overrep() {
  # args: fastqc_data.txt  -> prints "Sequence<TAB>Percentage"
  awk -v minpct="$MIN_PCT" '
    BEGIN{inmod=0}
    /^>>Overrepresented sequences/ {inmod=1; getline; next}
    /^>>END_MODULE/ {inmod=0}
    inmod && NR>1 {
      # columns: Sequence \t Count \t Percentage \t Possible Source
      seq=$1; pct=$3+0
      if (seq != "Sequence" && pct >= minpct) print seq "\t" pct
    }' "$1"
}

# ==== MAIN LOOP ====
shopt -s nullglob
for fq in *.fastq *.fastq.gz; do
  [ -e "$fq" ] || continue
  base="$(basename "$fq")"
  stem="${base%.fastq.gz}"
  stem="${stem%.fastq}"

  echo "==[1/5] FastQC: $base =="
  fastqc -t "$THREADS" -o qc --extract "$fq" > "logs/${stem}.fastqc.log" 2>&1

  # Find the extracted fastqc_data.txt (FastQC makes <stem>_fastqc/fastqc_data.txt)
  FQDIR_CAND1="qc/${stem}_fastqc/fastqc_data.txt"
  FQDIR_CAND2="qc/${base%.*}_fastqc/fastqc_data.txt"
  if   [ -f "$FQDIR_CAND1" ]; then FQDATA="$FQDIR_CAND1"
  elif [ -f "$FQDIR_CAND2" ]; then FQDATA="$FQDIR_CAND2"
  else
    echo "ERROR: fastqc_data.txt not found for $base" >&2
    continue
  fi

  echo "==[2/5] Parsing overrepresented sequences from FastQC =="
  mapfile -t overrep_lines < <(extract_overrep "$FQDATA")

  # Build per-sample adapters FASTA
  ADP="adapters/${stem}.fa"
  : > "$ADP"  # truncate
  # Always include the standard Illumina universal (common in library preps)
  echo -e ">illumina_universal\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCA" >> "$ADP"

  if $TRIM_HOMOPOLYMERS; then
    # Add homopolymer tails (helpful for poly-G/A artifacts)
    printf ">poly_G_%d\n" "$HOMOPOLYMER_LEN" >> "$ADP"
    printf "%0.sG" $(seq 1 "$HOMOPOLYMER_LEN") >> "$ADP"; echo >> "$ADP"
    printf ">poly_A_%d\n" "$HOMOPOLYMER_LEN" >> "$ADP"
    printf "%0.sA" $(seq 1 "$HOMOPOLYMER_LEN") >> "$ADP"; echo >> "$ADP"
  fi

  n_custom=0
  for line in "${overrep_lines[@]}"; do
    seq="${line%%$'\t'*}"
    # Skip if sequence is too short to be meaningful
    if [ "${#seq}" -lt 15 ]; then continue; fi
    n_custom=$((n_custom+1))
    printf ">overrep_%02d\n%s\n" "$n_custom" "$seq" >> "$ADP"
  done

  echo "   Added $n_custom sample-specific overrepresented sequences to $ADP"

  echo "==[3/5] Trimming with cutadapt =="
  # Build -a arguments from the adapters.fa
  # (cutadapt accepts multiple -a entries; we pass all)
  mapfile -t ADAPTERS < <(awk '/^>/ {next} {print}' "$ADP")
  args=()
  for a in "${ADAPTERS[@]}"; do
    args+=( -a "$a" )
  done

  out="trimmed/${stem}.overrep_trim.fastq.gz"
  cutadapt \
    -j "$THREADS" \
    "${args[@]}" \
    --trim-n \
    -m "$MIN_LEN" \
    -o "$out" \
    "$fq" \
    > "logs/${stem}.cutadapt.overrep.txt" 2>&1

  echo "==[4/5] Re-run FastQC on trimmed output =="
  fastqc -t "$THREADS" -o qc "$out" > "logs/${stem}.fastqc_posttrim.log" 2>&1

  echo "==[5/5] Done: $out =="
done

echo "All samples processed. Reports in ./logs and ./qc"
