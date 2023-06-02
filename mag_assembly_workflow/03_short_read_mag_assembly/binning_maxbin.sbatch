#!/bin/bash

#Slurm email messages
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-user=hibberdm@wustl.edu

ml metagenomics_pipeline/3.1
ml maxbin/2.2.5

source vars.sh

if [ -z "${1}" ]; then
    echo "Must pass the metadata"
    exit 1
fi

COUNTS_DIR="06_counts_subgroup"

BIN_DIR="09_contig_binning_maxbin"

mkdir -p $BIN_DIR

SUBGROUP=$( sed -ne '/^#/!p' "${1}" | sed -n "${SLURM_ARRAY_TASK_ID}"p )

mkdir -p ${BIN_DIR}/NZ_${SUBGROUP}

./depth_table_to_file.py -i contig_depths/NZ_${SUBGROUP}.depth -p contig_depths/

perl /opt/apps/labs/jglab/software/maxbin/2.2.7/run_MaxBin.pl \
    -contig ${ANNO_DIR}/${SUBGROUP}/NZ_${SUBGROUP}.fna \
    -out ${BIN_DIR}/NZ_${SUBGROUP}/NZ_${SUBGROUP}_bin \
    -abund_list contig_depths/NZ_${SUBGROUP}.depth_list \
    -thread ${SLURM_CPUS_PER_TASK} \
    -min_contig_length 1000 \
    -plotmarker \
    -markerset 107

