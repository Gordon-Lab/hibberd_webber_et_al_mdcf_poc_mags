#!/bin/bash

################################################################################
# SCRIPT NAME : megahit_array.sbatch                                           #
# DESCRIPTION : this script performs assembly using the megahit assembler      #
# ARGS        : mapping file (positional, optional)                            #
# AUTHOR      : Matthew C. Hibberd                                             #
# EMAIL       : hibberdm@wustl.edu                                             #
################################################################################

#Slurm email messages
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -e

if [ -z "${1}" ]; then
    echo "Must pass the metadata"
    exit 1
fi

ml metagenomics_pipeline/3.0

source vars.sh

mkdir -p "$ASSEMBLY_DIR"

SID=$( sed -ne '/^#/!p' "${1}" | sed -n "${SLURM_ARRAY_TASK_ID}"p )

echo -e "SID: ${SID}"

ml megahit/1.1.4

megahit \
    --num-cpu-threads "${SLURM_CPUS_PER_TASK}" \
    --memory 48000000000 \
    --mem-flag 1 \
    --verbose \
    -1 "${HF_DIR}"/concatenated/"${SID}"_HF_CAT_R1.fq.gz \
    -2 "${HF_DIR}"/concatenated/"${SID}"_HF_CAT_R2.fq.gz \
    --out-dir "${ASSEMBLY_DIR}/${SID}" \
    --out-prefix "$SID"
