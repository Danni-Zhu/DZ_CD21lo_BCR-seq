#!/bin/bash

for fq in /n/scratch3/users/y/yz348/CD21lo_BCRseq/results/merged_reads/QC_filtered/*.fastq

do


sbatch -p short -t 0-12:00 -c 20 --job-name PrimerMask_Dedup --mem 20G --wrap="sh /n/scratch3/users/y/yz348/CD21lo_BCRseq/scripts/Primer-mask_deDup.sh $fq"

sleep 1 # wait 1 second between each job submission


done
