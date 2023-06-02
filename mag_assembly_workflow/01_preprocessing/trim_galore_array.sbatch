#!/bin/bash

################################################################################
# SCRIPT NAME : trim_galore_array.sbatch                                       #
# DESCRIPTION : this script performs quality and adapter trimming              #
# ARGS        : mapping file (positional, required)                            #
# AUTHOR      : Matthew C. Hibberd                                             #
# EMAIL       : hibberdm@wustl.edu                                             #
################################################################################

#SBATCH --mem=1G
#SBATCH --cpus-per-task=4

set -e

if [ -z "${1}" ]; then
    echo "Must pass the mapping file"
    exit 1
fi

ml metagenomics_pipeline/3.0

source vars.sh

METADATA=$( sed -ne '/^#/!p' "${1}" | sed -n "${SLURM_ARRAY_TASK_ID}"p )

IFS=$'\t' read -r RUN_DIR ID SID SUBGROUP PLATFORM HOST <<< "$METADATA"
unset IFS

mkdir -p  "${QC_DIR}/${RUN_DIR}"

ml trim_galore/0.6.4

counter=0

for ((i = 0; i < 4; ++i)); do
    LANE=$(( i + 1 ))

    if [[ -f "${READS_DIR}/${RUN_DIR}/s_${LANE}_1_withindex_sequence.txt_${ID}.fq.gz" &&
    -f "${READS_DIR}/${RUN_DIR}/s_${LANE}_4_withindex_sequence.txt_${ID}.fq.gz" ]]; then

        trim_galore \
            --cores "${SLURM_CPUS_PER_TASK}" \
            --phred33 \
            --length 50 \
            -e 0.1 \
            --2colour 20 \
            --paired \
            --basename "${SID}_${LANE}" \
            --output_dir "${QC_DIR}/${RUN_DIR}" \
            --trim-n \
            --retain_unpaired \
            --length_1 75 \
            --length_2 75 \
            "${READS_DIR}/${RUN_DIR}/s_${LANE}_1_withindex_sequence.txt_${ID}.fq.gz" \
            "${READS_DIR}/${RUN_DIR}/s_${LANE}_4_withindex_sequence.txt_${ID}.fq.gz"

        counter=$((counter+1))

    else
        echo "${READS_DIR}/${RUN_DIR}/s_${LANE}_1_withindex_sequence.txt_${ID}.fq.gz or \
             ${READS_DIR}/${RUN_DIR}/s_${LANE}_4_withindex_sequence.txt_${ID}.fq.gz not found." 

    fi

done

if ! ((counter > 0)); then
    echo "[ERROR] No fastq.gz files found for ID ${ID} in ${READS_DIR}/${RUN_DIR}!"
    exit 1
fi

