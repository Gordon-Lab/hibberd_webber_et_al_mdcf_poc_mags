#!/bin/bash

#Slurm email messages
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

ml metagenomics_pipeline/3.1
ml metabat/2.12.1

source vars.sh

if [ -z "${1}" ]; then
    echo "Must pass the metadata"
    exit 1
fi

COUNTS_DIR="06_counts_subgroup"

BIN_DIR="09_contig_binning"

mkdir -p $BIN_DIR

SUBGROUP=$( sed -ne '/^#/!p' "${1}" | sed -n "${SLURM_ARRAY_TASK_ID}"p )

mkdir -p ${BIN_DIR}/NZ_${SUBGROUP}

jgi_summarize_bam_contig_depths \
    --outputDepth ${BIN_DIR}/NZ_${SUBGROUP}/NZ_${SUBGROUP}.depth \
    --pairedContigs ${BIN_DIR}/NZ_${SUBGROUP}/NZ_${SUBGROUP}.paired \
    --minContigLength 1000 \
    --minContigDepth 1 \
    ${COUNTS_DIR}/NZ_${SUBGROUP}*_nodup.bam

metabat2 \
    --unbinned \
    --saveCls \
    --inFile ${ANNO_DIR}/${SUBGROUP}/NZ_${SUBGROUP}.fna \
    --outFile ${BIN_DIR}/NZ_${SUBGROUP}/NZ_${SUBGROUP}_bin \
    --abdFile ${BIN_DIR}/NZ_${SUBGROUP}/NZ_${SUBGROUP}.depth \
    -v

