#!/bin/bash

#Slurm email messages
#SBATCH --mem=36G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-user="hibberdm@wustl.edu"
#SBATCH -o slurm-checkm-%A.out

set -e

# spack load CheckM 1.1.3
eval $(spack load --sh py-checkm@1.1.3)

# run checkm to extract and align marker gene sequences
checkm \
    lineage_wf \
    genomes \
    drep_output \
    -t $SLURM_CPUS_PER_TASK












 
  
   
    
