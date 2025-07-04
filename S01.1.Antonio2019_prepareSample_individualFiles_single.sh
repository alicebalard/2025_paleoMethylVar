#!/bin/bash
#$ -N prepAntonio
#$ -S /bin/bash
#$ -l h_vmem=12G,tmem=12G ## NB: tmem value is per core
#$ -t 1-134
#$ -tc 100 # 100 samples max in one go
#$ -pe smp 4 # Request N cores per task 
#$ -l h_rt=72:00:00
#$ -wd /SAN/ghlab/epigen/Alice/paleo_project/logs # one err and out file per sample
#$ -R y # reserve the resources, i.e. stop smaller jobs from getting into the queue while you wait for all the required resources to become available for you
THREADS=4

##########################################
## Download hs37d5 reference and index it:

## hs37d5: Includes data from GRCh37, the rCRS mitochondrial sequence, Human herpesvirus 4 type 1 and the concatenated decoy sequences. Data is in one file, hs37d5.fa.gz, hosted by the EBI FTP site
## wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
## gunzip reference/hs37d5.fa.gz ## done

# Requirements: bwa, samtools, pmdtools
BWA="/share/apps/genomics/bwa-0.7.19/bwa"
# source python 3.6.4 to run cutadapt 2.4
source /share/apps/source_files/python/python-3.6.4.source
# add FastQC to my path
export PATH=/share/apps/genomics/FastQC-0.11.9/:$PATH
# add pigz to my path
export PATH=/share/apps/pigz-2.6/:$PATH
## $BWA index reference/hs37d5.fa ## done

DATADIR="/SAN/ghlab/epigen/Alice/paleo_project/data"
SAMPLEDIR="/SAN/ghlab/epigen/Alice/paleo_project/data/01rawfastq/Antonio2019"
REFERENCE="$DATADIR/reference/hs37d5.fa"

## We make a temporary folder for all intermediate files that we delete at the end.
TEMP_OUTDIR="$SAMPLEDIR/TEMP"

mkdir -p $TEMP_OUTDIR
cd $TEMP_OUTDIR

# Create the files to loop over if it does not exist:
ls -1 $SAMPLEDIR/*fastq.gz > $TEMP_OUTDIR/list_of_files.tmp
xargs -n1 basename < $TEMP_OUTDIR/list_of_files.tmp | cut -d. -f1 > $TEMP_OUTDIR/list_of_sample_names.tmp

## Select the correct line of list of files at each iteration
FASTQ=$(sed -n "${SGE_TASK_ID}p" $TEMP_OUTDIR/list_of_files.tmp)
SAMPLE_NAME=$(sed -n "${SGE_TASK_ID}p" $TEMP_OUTDIR/list_of_sample_names.tmp)

# Log the start of the job. 
echo "**** Job $JOB_NAME.$JOB_ID started at $(date) ****"
echo "Task ID: $SGE_TASK_ID"
echo "Selected input file: $SAMPLE_NAME"

TRIMOM_OUTDIR="$TEMP_OUTDIR/01Trimmed_data"
mkdir -p $TRIMOM_OUTDIR

if [ ! -f "$DATADIR/02samplesBams/Antonio2019/${SAMPLE_NAME}.filtered.sorted.dedup.bam" ]; then
    
    ## Trim raw reads
    ### Single reads: "Raw reads were trimmed using Trim Galore with no quality filter and a length filter of 30 bp (-q0, --length 30, -a ‘AGATCGGAAGAGCACACGTCTGAACTCC’)."

    TRIMMED="$TRIMOM_OUTDIR/${SAMPLE_NAME}_trimmed.fq.gz"

    echo "**** Start of step 1: Trim galore with with no quality filter and a length filter of 30 bp $(date) ****" 
    echo "Run Trim Galore command without fastqc"
    /share/apps/genomics/TrimGalore-0.6.7/trim_galore --path_to_cutadapt /share/apps/genomics/cutadapt-2.5/bin/cutadapt -q 0 --length 30 -a AGATCGGAAGAGCACACGTCTGAACTCC --cores 1 --output_dir $TRIMOM_OUTDIR $FASTQ ## the help advise for 4 cores, but Ed Martin from CS advise that with pigz, 1 core is better on this system. Or we could use bgzip

    ## Pooja:
    #for fastq in $(ls *fastq.gz)
    #do
    #name=$(basename $fastq | sed 's/.fastq.gz//g')
    #AdapterRemoval \
	#--file1 $fastq \
	#--trimns --trimqualities --gzip --minlength 35 --minquality 20 \
	#--adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
	#--adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	#--basename $name
    #done

    echo "**** End of step 1: $(date) ****"

    ################################
    # Ancient DNA-optimized Mapping:
    ## map reads to Homo sapiens genome assembly hs37d5

    OUTPUT_SAM="$TRIMOM_OUTDIR/${SAMPLE_NAME}.sam"

    echo "Align sample ${INPUT##*/} to hs37d5..." ## 2 days to run
    $BWA aln -l 1024 -n 0.01 -o 2 -t $THREADS "$REFERENCE" "$TRIMMED" > "$TRIMMED.sai" ## (long seed, high mismatch tolerance)
    ## Updated after to see if difference "Improving ancient DNA read mapping against modern reference genomes - BMC Genomics"
    ## Could improve methylation

    echo "Generate alignment file (single ended):"
    $BWA samse "$REFERENCE" "$TRIMMED.sai" "$TRIMMED" > "$OUTPUT_SAM"

    ## Çokoğlu et al 2024 "We filtered out reads of size less than 35 bps, with a mapping quality (MAPQ) of less than 30, and with more than 10% mismatches to the reference genome."
    # 1. Filter by MAPQ
    samtools view -b -q 30 -o "$TRIMOM_OUTDIR/${SAMPLE_NAME}.mapq30.bam" "$OUTPUT_SAM"

    # 2. Filter by length and mismatch fraction (retain header)
    samtools view -h "$TRIMOM_OUTDIR/${SAMPLE_NAME}.mapq30.bam" | \
	awk 'BEGIN {OFS="\t"} {
  if ($0 ~ /^@/) { print; next }  # keep header
  seq_len = length($10);
  nm = -1;
  for(i=12; i<=NF; i++) {
    if($i ~ /^NM:i:/) {
      split($i, a, ":");
      nm = a[3];
      break;
    }
  }
  if (nm == -1) next;  # skip if NM tag not found
  mismatch_frac = nm / seq_len;
  if(seq_len >= 35 && mismatch_frac <= 0.10) print;
}' | samtools view -b -o "$TRIMOM_OUTDIR/${SAMPLE_NAME}.filtered.bam"

    # 3. Sort and index
    samtools sort -o "$TRIMOM_OUTDIR/${SAMPLE_NAME}.filtered.sorted.bam" "$TRIMOM_OUTDIR/${SAMPLE_NAME}.filtered.bam"
    samtools index "$TRIMOM_OUTDIR/${SAMPLE_NAME}.filtered.sorted.bam"

    echo "Duplicates removal:"
    java="/share/apps/java/bin/java"
    "$java" -jar -Xmx3G -Xms3G /share/apps/genomics/picard-2.20.3/bin/picard.jar MarkDuplicates \
	    I="$TRIMOM_OUTDIR/${SAMPLE_NAME}.filtered.sorted.bam" \
	    O="$TRIMOM_OUTDIR/${SAMPLE_NAME}.filtered.sorted.dedup.bam" \
	    VALIDATION_STRINGENCY=LENIENT \
	    REMOVE_DUPLICATES=TRUE \
	    M="$TRIMOM_OUTDIR/${SAMPLE_NAME}.filtered.sorted.bam.marked_dup_metrics.txt"

    echo "Move final file in the correct directory:"
    if [ ! -f "$DATADIR/02samplesBams/Antonio2019/${SAMPLE_NAME}.filtered.sorted.dedup.bam" ]; then
	mv "$TRIMOM_OUTDIR/${SAMPLE_NAME}.filtered.sorted.dedup.bam" "$DATADIR/02samplesBams/Antonio2019/"
    fi   
    echo "Remove the temporary files:"
    rm "$TRIMOM_OUTDIR/${SAMPLE_NAME}."*
    rm "$TRIMOM_OUTDIR/${SAMPLE_NAME}_"*
else
    echo "Sample already processed"
    echo "Move final file in the correct directory if needed:"
    if [ ! -f "$DATADIR/02samplesBams/Antonio2019/${SAMPLE_NAME}.filtered.sorted.dedup.bam" ]; then
	mv "$TRIMOM_OUTDIR/${SAMPLE_NAME}.filtered.sorted.dedup.bam" "$DATADIR/02samplesBams/Antonio2019/"
    fi
    echo "Remove the temporary files:"
    rm "$TRIMOM_OUTDIR/${SAMPLE_NAME}."*
    rm "$TRIMOM_OUTDIR/${SAMPLE_NAME}_"*
fi    
echo "Analysis done!"
