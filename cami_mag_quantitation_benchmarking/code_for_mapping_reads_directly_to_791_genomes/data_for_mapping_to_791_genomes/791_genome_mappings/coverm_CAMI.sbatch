#!/bin/bash
#SBATCH --mem=50G
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-user="hmbucklin@wustl.edu"
#SBATCH --array=0-63

export PATH=$PATH:/scratch/jglab/hmbucklin/CoverM_download2/coverm-x86_64-unknown-linux-musl-0.6.1.tar.gz

coverm genome \
--bam-files /scratch/jglab/hmbucklin/03032021_CAMI_mouse_reanalysis/redo_corrected_bowtie_counts_all_samples/sample_"${SLURM_ARRAY_TASK_ID}"_nodup.bam \
--genome-fasta-directory /scratch/jglab/hmbucklin/03032021_CAMI_mouse_reanalysis/791_genomes_fna \
--methods rpkm \
> sample_"${SLURM_ARRAY_TASK_ID}"_coverage_rpkm.txt
