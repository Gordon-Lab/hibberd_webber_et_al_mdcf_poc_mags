#!/bin/bash

################################################################################
# SCRIPT NAME : cat_subgroup_array.sbatch                                      #
# DESCRIPTION : this script merges data between sequencing runs by SID         #
# ARGS        : mapping file (positional, required)                            #
# AUTHOR      : Matthew C. Hibberd                                             #
# EMAIL       : hibberdm@wustl.edu                                             #
################################################################################

set -e

if [ -z "${1}" ]; then
    echo "Must pass the metadata"
    exit 1
fi

ml metagenomics_pipeline/3.0

source vars.sh

mkdir -p "${HF_DIR}/subgroup"

SID=$( sed -ne '/^#/!p' "${1}" | sed -n "${SLURM_ARRAY_TASK_ID}"p )

echo "${SID}"

# R1 File Concatenation

R1_FILES=$(ls ${HF_DIR}/*/${SID}_*_HF.fq.1.gz)
echo -e "R1_FILES:\n${R1_FILES}"

cat \
    ${HF_DIR}/*/${SID}_*_HF.fq.1.gz \
    > "${HF_DIR}"/concatenated/"${SID}"_HF_CAT_R1.fq.gz

# R2 File Concatenation

R2_FILES=$(ls ${HF_DIR}/*/${SID}_*_HF.fq.2.gz)
echo -e "R2_FILES:\n${R2_FILES}"

cat \
    ${HF_DIR}/*/${SID}_*_HF.fq.2.gz \
    > "${HF_DIR}"/concatenated/"${SID}"_HF_CAT_R2.fq.gz


#rm ${R1_FILES}
#rm ${R2_FILES}

