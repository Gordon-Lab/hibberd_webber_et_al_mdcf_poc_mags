#!/bin/bash
#SBATCH --cpus-per-task=24
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
#SBATCH --mail-user="hmbucklin@wustl.edu"
#SBATCH --mem=100G

ml kallisto/0.43.0

kallisto index \
-i kalliso_791_Ns_CAMI_index \
/scratch/jglab/hmbucklin/03032021_CAMI_mouse_reanalysis/791_genomes_w_Ns/*.fna