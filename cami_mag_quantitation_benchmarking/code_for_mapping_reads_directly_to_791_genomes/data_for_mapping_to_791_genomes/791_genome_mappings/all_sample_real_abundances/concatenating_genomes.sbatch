#!/bin/bash
#SBATCH --mem=10G
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
#SBATCH --mail-user="hmbucklin@wustl.edu"

cat /scratch/jglab/hmbucklin/03032021_CAMI_mouse_reanalysis/791_genomes/*.fa > all_791_genomes_concatenated.fa