#!/bin/bash
#SBATCH --mem=10G
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
#SBATCH --mail-user="hmbucklin@wustl.edu"
#SBATCH --array=0-63

gunzip -c /scratch/jglab/hmbucklin/02092021_CAMI_mouse_data_resync/19122017_mousegut_scaffolds/2017.12.29_11.37.26_sample_"${SLURM_ARRAY_TASK_ID}"/reads/anonymous_reads.fq.gz > /scratch/jglab/hmbucklin/02092021_CAMI_mouse_data_resync/19122017_mousegut_scaffolds/2017.12.29_11.37.26_sample_"${SLURM_ARRAY_TASK_ID}"/reads/sample_"${SLURM_ARRAY_TASK_ID}"_all_reads.fq

paste - - - - - - - - < /scratch/jglab/hmbucklin/02092021_CAMI_mouse_data_resync/19122017_mousegut_scaffolds/2017.12.29_11.37.26_sample_"${SLURM_ARRAY_TASK_ID}"/reads/sample_"${SLURM_ARRAY_TASK_ID}"_all_reads.fq \
    | tee >(cut -f 1-4 | tr "\t" "\n" > /scratch/jglab/hmbucklin/02092021_CAMI_mouse_data_resync/19122017_mousegut_scaffolds/2017.12.29_11.37.26_sample_"${SLURM_ARRAY_TASK_ID}"/reads/sample_"${SLURM_ARRAY_TASK_ID}"_R1_reads.fq) \
    |       cut -f 5-8 | tr "\t" "\n" > /scratch/jglab/hmbucklin/02092021_CAMI_mouse_data_resync/19122017_mousegut_scaffolds/2017.12.29_11.37.26_sample_"${SLURM_ARRAY_TASK_ID}"/reads/sample_"${SLURM_ARRAY_TASK_ID}"_R2_reads.fq