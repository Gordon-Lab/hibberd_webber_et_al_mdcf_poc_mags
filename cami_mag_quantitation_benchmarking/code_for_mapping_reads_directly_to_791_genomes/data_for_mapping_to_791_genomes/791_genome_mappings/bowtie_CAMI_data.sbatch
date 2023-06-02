#!/bin/bash
#SBATCH --mem=100G
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-user="hmbucklin@wustl.edu"
#SBATCH --array=0-63

#Reading in modules needed for the script
ml bowtie2/2.3.4.1
ml samtools/1.9
ml bedtools/2.27.1
ml subread/1.6.2

COUNTS_DIR="/scratch/jglab/hmbucklin/03032021_CAMI_mouse_reanalysis/bowtie_counts_all_samples"
SAMPLE_NAME=sample_"${SLURM_ARRAY_TASK_ID}"

if [ -f /scratch/jglab/hmbucklin/03032021_CAMI_mouse_reanalysis/791_genomes/all_791_genomes_concatenated.fa ]; then

        if [ ! -f /scratch/jglab/hmbucklin/03032021_CAMI_mouse_reanalysis/791_genomes/all_791_genomes_concatenated.fa.1.bt2 ]; then
        bowtie2-build /scratch/jglab/hmbucklin/03032021_CAMI_mouse_reanalysis/791_genomes/all_791_genomes_concatenated.fa /scratch/jglab/hmbucklin/03032021_CAMI_mouse_reanalysis/791_genomes/all_791_genomes_concatenated.fa
    fi
fi

   bowtie2 \
        -p "${SLURM_CPUS_PER_TASK}" \
        -1 /scratch/jglab/hmbucklin/02092021_CAMI_mouse_data_resync/19122017_mousegut_scaffolds/2017.12.29_11.37.26_sample_"${SLURM_ARRAY_TASK_ID}"/reads/sample_"${SLURM_ARRAY_TASK_ID}"_R1_reads.fq \
        -2 /scratch/jglab/hmbucklin/02092021_CAMI_mouse_data_resync/19122017_mousegut_scaffolds/2017.12.29_11.37.26_sample_"${SLURM_ARRAY_TASK_ID}"/reads/sample_"${SLURM_ARRAY_TASK_ID}"_R2_reads.fq \
        -x /scratch/jglab/hmbucklin/03032021_CAMI_mouse_reanalysis/791_genomes/all_791_genomes_concatenated.fa \
        | samtools view --threads "${SLURM_CPUS_PER_TASK}" -huS - \
        | samtools sort --threads "${SLURM_CPUS_PER_TASK}" -n - \
        | samtools fixmate --threads "${SLURM_CPUS_PER_TASK}" -m - - \
        | samtools sort --threads "${SLURM_CPUS_PER_TASK}" - \
        | samtools markdup --threads "${SLURM_CPUS_PER_TASK}" -rs - "${COUNTS_DIR}/${SAMPLE_NAME}_nodup.bam"

         #header would be "chromosome" "depth" "bases with depth > 2" "chromosome size" "fraction of bases on chromosome with depth = depth"        

