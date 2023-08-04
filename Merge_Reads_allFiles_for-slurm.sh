#!/bin/bash

for fq in /n/scratch3/users/y/yz348/CD21lo_BCRseq/raw_data/FC_07903/Unaligned_1a_PF_mm0/Data/Project_dannizhu/*R1*.fastq

do


sbatch -p short -t 0-6:00 -c 12 --job-name Merge_pairedEnd_Reads --mem 100G --wrap="sh /n/scratch3/users/y/yz348/CD21lo_BCRseq/scripts/Merge_Reads.sh $fq"

sleep 1 # wait 1 second between each job submission


done
