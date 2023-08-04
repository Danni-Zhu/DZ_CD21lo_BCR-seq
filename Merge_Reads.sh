#!/bin/bash

# This script takes a fastq file with paired-end reads of BCR-seq data and assembles paired-end reads into a complete sequence with pRESTO 

# initialize variable to store name of input fastq file 
fq=$1

# grab base of filename for naming outputs 

samplename=`basename $fq R1_001.fastq`
echo "Sample name is $samplename" 
read_one_fq=${samplename}R1_001.fastq
read_two_fq=${samplename}R2_001.fastq

path_to_fq=/n/scratch3/users/y/yz348/CD21lo_BCRseq/raw_data/FC_07903/Unaligned_1a_PF_mm0/Data/Project_dannizhu/
merge_out=/n/scratch3/users/y/yz348/CD21lo_BCRseq/results/merged_reads/${samplename}_
path_out=/n/scratch3/users/y/yz348/CD21lo_BCRseq/results/merged_reads/

echo "Processing file %fq"

module load gcc/9.2.0 python/3.9.14

#run pRESTO

AssemblePairs.py align -1 $path_to_fq$read_one_fq -2 $path_to_fq$read_two_fq \
--coord illumina \
--rc tail \
--outname $merge_out \
--outdir $path_out




