cd /mnt/c/Users/hamda/Desktop/UGI_thingy/3rdyear/genome/quality

for f in *.fastq; do
  base=${f%.fastq}
  echo "Processing $f ..."
  cutadapt \
    -j 8 \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \  # !!!!!! 
    -m 30 \
    --trim-n \
    -o trimmed/${base}.trim.fastq.gz \
    "$f" > logs/${base}.cutadapt.txt
done


