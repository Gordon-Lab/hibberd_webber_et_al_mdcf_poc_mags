#!/bin/bash

#Slurm email messages
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-user="hibberdm@wustl.edu"
#SBATCH --mem=16G
#SBATCH --cpus-per-task=16

ml py-drep/2.3.2-python-3.6.5
ml prodigal/2.6.3
ml mash/2.2
ml checkm/1.1.3

dRep dereplicate \
  drep_output_99_v1 \
  -g final/*.mod2.fa \
  -p ${SLURM_CPUS_PER_TASK} \
  -l 50000 \
  --completeness 90 \
  --contamination 5 \
  --S_algorithm ANImf \
  --P_ani 0.9 \
  --S_ani 0.99 \
  --cov_thresh 0.1 \
  --genomeInfo 210608_MDCF_POC_MAG_checkm.csv

