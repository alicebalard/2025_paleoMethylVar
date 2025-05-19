#!/bin/bash
#$ -N copyfiles
#$ -S /bin/bash
#$ -l tmem=5G
#$ -l h_vmem=5G
#$ -l h_rt=24:00:00
#$ -wd /SAN/ghlab/epigen/Alice/paleo_project/logs # one err and out file per sample
#$ -R y # reserve the resources, i.e. stop smaller jobs from getting into the queue while you wait for all the required resources to become available for you

