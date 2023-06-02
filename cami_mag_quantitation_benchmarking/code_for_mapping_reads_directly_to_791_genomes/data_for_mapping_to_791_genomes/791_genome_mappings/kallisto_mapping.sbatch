#!/bin/bash
#SBATCH --cpus-per-task=24
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-user="hmbucklin@wustl.edu"
#SBATCH --mem=100G
#SBATCH --array=0-63

ml kallisto/0.43.0

kallisto quant \
-i kalliso_791_Ns_CAMI_index \
-b 100 \
-o sample_"${SLURM_ARRAY_TASK_ID}" \
/scratch/jglab/hmbucklin/02092021_CAMI_mouse_data_resync/19122017_mousegut_scaffolds/2017.12.29_11.37.26_sample_"${SLURM_ARRAY_TASK_ID}"/reads/sample_"${SLURM_ARRAY_TASK_ID}"_R1_reads.fq \
/scratch/jglab/hmbucklin/02092021_CAMI_mouse_data_resync/19122017_mousegut_scaffolds/2017.12.29_11.37.26_sample_"${SLURM_ARRAY_TASK_ID}"/reads/sample_"${SLURM_ARRAY_TASK_ID}"_R2_reads.fq
