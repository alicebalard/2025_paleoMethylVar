#!/bin/bash
#$ -N paleo_generateBedHg19
#$ -S /bin/bash
#$ -pe smp 4 # Request N cores per task                                                                                                                                                      
#$ -l tmem=8G
#$ -l h_vmem=8G
#$ -l h_rt=48:00:00
#$ -wd /SAN/ghlab/pophistory/Alice/paleo_project/logs # one err and out file per sample
#$ -R y # reserve the resources, i.e. stop smaller jobs from getting into the queue while you wait for all the required resources to become available for you

Rscript /SAN/ghlab/pophistory/Alice/paleo_project/code/generatehg19CpGs.R
