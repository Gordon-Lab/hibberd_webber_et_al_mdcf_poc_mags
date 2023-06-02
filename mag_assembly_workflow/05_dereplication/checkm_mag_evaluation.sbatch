#!/bin/bash

#Slurm email messages
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4

set -e

ml hmmer/3.1b2
ml pplacer/1.1.alpha19
ml prodigal/2.6.2

export PATH=$PATH:/scratch/jglab/hibberdm/checkm/bin/
export PYTHONPATH=$PYTHONPATH:/scratch/jglab/hibberdm/checkm/

#checkm data setRoot /scratch/ref/checkm/

#source vars.sh

if [ -z "${1}" ]; then
    echo "Must pass the metadata"
    exit 1
fi

ANNO_DIR="bins"

MAG=$( sed -ne '/^#/!p' "${1}" | sed -n "${SLURM_ARRAY_TASK_ID}"p )

checkm \
    lineage_wf \
    ${ANNO_DIR} \
    ${ANNO_DIR}/${MAG}_checkm_output \
    -x ${MAG}.fa \
    -f ${ANNO_DIR}/${MAG}_checkm_results.txt \
    -t ${SLURM_CPUS_PER_TASK} \
    --tab_table

#ls ${COUNTS_DIR}/NZ_${SUBGROUP}*_nodup.bam | xargs -n1 -P5 samtools index

#checkm \
#    coverage \
#    ${BIN_DIR}/NZ_${SUBGROUP} \
#    ${BIN_DIR}/NZ_${SUBGROUP}/NZ_${SUBGROUP}.coverage \
#    ${COUNTS_DIR}/NZ_${SUBGROUP}*_nodup.bam \
#    -x fa \
#    -t ${SLURM_CPUS_PER_TASK}

checkm \
    tetra \
    ${ANNO_DIR}/${MAG}.fa \
    ${ANNO_DIR}/${MAG}.tetra

checkm \
    outliers \
    ${ANNO_DIR}/${MAG}_checkm_output \
    ${ANNO_DIR}/ \
    ${ANNO_DIR}/${MAG}.tetra \
    ${ANNO_DIR}/${MAG}.outliers \
    -x ${MAG}.fa

checkm \
    modify \
    ${ANNO_DIR}/${MAG}.fa \
    ${ANNO_DIR}/${MAG}.fa \
    ${ANNO_DIR}/${MAG}.mod.fa \
    --outlier_file ${ANNO_DIR}/${MAG}.outliers

checkm \
    lineage_wf \
    ${ANNO_DIR} \
    ${ANNO_DIR}/${MAG}_checkm_output2 \
    -x ${MAG}.mod.fa \
    -f ${ANNO_DIR}/${MAG}_checkm_results_mod.txt \
    -t ${SLURM_CPUS_PER_TASK} \
    --tab_table

#checkm \
#    profile \
#    ${BIN_DIR}/NZ_${SUBGROUP}/NZ_${SUBGROUP}.coverage \
#    -f ${BIN_DIR}/NZ_${SUBGROUP}/NZ_${SUBGROUP}.bin_profile \
#    --tab_table

#checkm \
#    ssu_finder \
#    ${ANNO_DIR}/${SUBGROUP}/NZ_${SUBGROUP}.fna \
#    ${BIN_DIR}/NZ_${SUBGROUP} \
#    ${BIN_DIR}/NZ_${SUBGROUP}/checkm_ssu_finder \
#    -x fa

#if [ ! -f "${BIN_DIR}/bacteria.ms" ]; then
#    checkm taxon_set domain Bacteria ${BIN_DIR}/bacteria.ms
#fi

#checkm \
#    merge \
#    ${BIN_DIR}/bacteria.ms \
#    ${BIN_DIR}/NZ_${SUBGROUP} \
#    ${BIN_DIR}/NZ_${SUBGROUP}/checkm_merge \
#    -x fa \
#    -t ${SLURM_CPUS_PER_TASK}

