#!/usr/bin/env bash

#SBATCH --mem 80G
#SBATCH --cpus-per-task 8
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-user="hibberdm@wustl.edu"

ml kallisto/0.43.0

export WORK_DIR="/scratch/jglab/hibberdm/metagenomic_assembly/200213_MDCF_POC"
export HF_DIR=${WORK_DIR}/02_host_filtered
COUNTS_DIR="mag_ab_pseudo_checkm_magpurify_clean_kraken2clean"
IDX="all_bins.kallisto"

mkdir -p "$COUNTS_DIR"

METADATA=$( sed -ne '/^#/!p' "${1}" | sed -n "${SLURM_ARRAY_TASK_ID}"p )

IFS=$'\t' read -r SID SUBGROUP <<< "$METADATA"
unset IFS

prefix="NZ_"${SID}

kallisto \
  quant \
  -b 100 \
  -i $IDX \
  -o ${COUNTS_DIR}/$prefix \
  -t $SLURM_CPUS_PER_TASK \
  "${HF_DIR}"/concatenated/"${SID}"_HF_CAT_R1.fq.gz \
  "${HF_DIR}"/concatenated/"${SID}"_HF_CAT_R2.fq.gz
