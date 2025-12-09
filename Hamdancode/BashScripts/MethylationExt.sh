cd /mnt/c/Users/hamda/Desktop/UGI_thingy/3rdyear/genome/quality

BAM="align_hg38/SP75_bisulfite.trim_bismark_bt2.deduplicated.bam"  # pick one that exists
bismark_methylation_extractor \
  --gzip --bedGraph --comprehensive \
  --parallel 8 --buffer_size 6G \
  --genome_folder /mnt/c/Users/hamda/Desktop/UGI_thingy/3rdyear/genome/bis/allgenomes/refs/hg38 \
  -o methylation_hg38 \
  "$BAM" 2>&1 | tee methylation_hg38/SP75.extract.log

ls -lh methylation_hg38 | head


cd /mnt/c/Users/hamda/Desktop/UGI_thingy/3rdyear/genome/quality

bismark_methylation_extractor \
  --gzip --bedGraph --comprehensive \
  --parallel 8 --buffer_size 6G \
  --genome_folder /mnt/c/Users/hamda/Desktop/UGI_thingy/3rdyear/genome/bis/allgenomes/refs/hg38 \
  -o methylation_hg38 \
  align_hg38/Zvej16_bisulfite.trim_bismark_bt2.deduplicated.bam \
  2>&1 | tee methylation_hg38/Zvej16.extract.log
