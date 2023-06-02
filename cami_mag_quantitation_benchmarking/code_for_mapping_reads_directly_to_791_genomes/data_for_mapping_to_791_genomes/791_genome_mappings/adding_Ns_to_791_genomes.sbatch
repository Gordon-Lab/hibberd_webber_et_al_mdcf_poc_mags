#!/bin/bash
#SBATCH --cpus-per-task=24
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
#SBATCH --mail-user="hmbucklin@wustl.edu"
#SBATCH --mem=50G

ml py-biopython/1.73-python-3.6.5

./kallisto_index_fasta_preprocess.py \
-d /scratch/jglab/hmbucklin/03032021_CAMI_mouse_reanalysis/791_genomes \
-o /scratch/jglab/hmbucklin/03032021_CAMI_mouse_reanalysis/791_genomes_w_Ns \
-e .fa
