#!/bin/bash

#Slurm email messages
#SBATCH --mem=2G
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-user="hibberdm@wustl.edu"
#SBATCH -o slurm-fasttree-%A.out

set -e

# load fasttree 2.1.10
eval $(spack load --sh fasttree@2.1.10)

# settings for multithreading
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# variable for version tracking
DATECODE=$( date +"%y%m%d_%H%M" )

# define inputs
FASTA="aligned_marker_gene_fasta_from_checkm.fasta"
TREE="${DATECODE}_output_tree.fasttree"

# run fasttree
FastTreeMP \
    -log "${DATECODE}_fasttree.log" \
    -gtr \
    -gamma \
     $FASTA \
    > $TREE
    
