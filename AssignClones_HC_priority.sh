#!/bin/bash

#SBATCH -c 20
#SBATCH -t 2-12:00
#SBATCH -p priority
#SBATCH -o DefineClone_%j.out
#SBATCH --job-name DefineClone
#SBATCH --mem 100G

DefineClones.py -d /n/scratch3/users/y/yz348/CD21lo_BCRseq/results/IMGT/HC_LC_separate/Reprocessed_HC_Combined_dist_ham.tsv --act set --model ham --norm len --dist 0.045


