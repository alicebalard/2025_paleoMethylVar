for fq in *.trim.fastq *.trim.fastq.gz; do
  [ -e "$fq" ] || continue
  echo "Mapping $fq ..."
  bismark \
    --genome "$HG38" \
    --non_directional \
    --score_min L,0,-0.6 \
    --bowtie2 -p 8 \
    --temp_dir tmp \
    -o align_hg38 \
    "$fq"
done
