#!/usr/bin/env bash

#SBATCH --mem 90G
#SBATCH --cpus-per-task 12
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-user="hibberdm@wustl.edu"

ml kallisto/0.43.0

cat collapsed/*.fna > all_bins.fna

kallisto \
  index \
  -i all_bins.kallisto \
  all_bins.fna
