#!/bin/bash

# This script takes in quality filtered, merged pair-end BCR-seq reads and remove the V-segment and C-region PCR primers. The script then collapse duplicated sequences and filter to unique sequences. The output files are ready for submission to IMGT.  

# initialize variable to store name of input fastq file 
fq=$1

# grab base of filename for naming outputs 

samplename=`basename $fq __quality-pass.fastq`
echo "Sample name is $samplename" 

# define variable for path to QC filtered fastq files 
PATH_TO_QCfastq=/n/scratch3/users/y/yz348/CD21lo_BCRseq/results/merged_reads/QC_filtered

mkdir -p $PATH_TO_QCfastq/primer_masked
mkdir -p $PATH_TO_QCfastq/deduplicated
mkdir -p $PATH_TO_QCfastq/log

### Masking primers 
MaskPrimers.py align -s $PATH_TO_QCfastq/${samplename}__quality-pass.fastq -p $PATH_TO_QCfastq/DZ_Vprimer.fasta \
--mode mask --pf VPRIMER --outname ${samplename}-FWD --log MPV.log --failed --outdir $PATH_TO_QCfastq/primer_masked
MaskPrimers.py align -s $PATH_TO_QCfastq/primer_masked/${samplename}-FWD_primers-pass.fastq -p $PATH_TO_QCfastq/DZ_Cprimer.fasta \
--mode cut --revpr --pf CPRIMER --outname ${samplename}-REV --log MPC.log --failed --outdir $PATH_TO_QCfastq/primer_masked
### create log file for primer masking 
ParseLog.py -l MPV.log MPC.log -f ID PRIMER ERROR --outdir $PATH_TO_QCfastq/log

### Collapse duplicated sequences 
CollapseSeq.py -s $PATH_TO_QCfastq/primer_masked/${samplename}-REV_primers-pass.fastq -n 20 --inner --uf CPRIMER \
    --cf VPRIMER --outdir $PATH_TO_QCfastq/deduplicated --act set --outname ${samplename}
### Filtering to repeated sequences, subset data to only unique sequences with at least two representative reads  
SplitSeq.py group -s $PATH_TO_QCfastq/deduplicated/${samplename}_collapse-unique.fastq -f DUPCOUNT --num 2 --outname ${samplename} --outdir $PATH_TO_QCfastq/deduplicated

### Create annotation table 
ParseHeaders.py table -s $PATH_TO_QCfastq/deduplicated/${samplename}_atleast-2.fastq -f ID DUPCOUNT CPRIMER VPRIMER --outdir $PATH_TO_QCfastq/deduplicated

