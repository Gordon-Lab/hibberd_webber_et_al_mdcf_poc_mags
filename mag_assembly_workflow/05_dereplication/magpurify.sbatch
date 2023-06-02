#!/bin/bash

#Slurm email messages
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4

set -e

ml hmmer/3.1b2
ml prodigal/2.6.2
ml ncbi-blast/2.6.0+
ml last/926
ml mash/2.2
# also coverm?
ml pplacer/1.1.alpha19

export PATH=$PATH:/scratch/jglab/hibberdm/checkm/bin/
export PYTHONPATH=$PYTHONPATH:/scratch/jglab/hibberdm/checkm/

export PATH=$PATH:/scratch/jglab/hibberdm/magpurify/bin/
export PYTHONPATH=$PYTHONPATH:/scratch/jglab/hibberdm/magpurify/
export MAGPURIFYDB=/scratch/jglab/hibberdm/magpurify/MAGpurify-db-v1.0

ANNO_DIR="bins"

MAG=$( sed -ne '/^#/!p' "${1}" | sed -n "${SLURM_ARRAY_TASK_ID}"p )

magpurify \
    phylo-markers \
    ${ANNO_DIR}/${MAG}.mod.fa \
    ${ANNO_DIR}/${MAG}_magpurify_output

magpurify \
    clade-markers \
    ${ANNO_DIR}/${MAG}.mod.fa \
    ${ANNO_DIR}/${MAG}_magpurify_output

magpurify \
    clean-bin \
    ${ANNO_DIR}/${MAG}.mod.fa \
    ${ANNO_DIR}/${MAG}_magpurify_output \
    ${ANNO_DIR}/${MAG}.mod2.fa

checkm \
    lineage_wf \
    ${ANNO_DIR} \
    ${ANNO_DIR}/${MAG}_checkm_output3 \
    -x ${MAG}.mod2.fa \
    -f ${ANNO_DIR}/${MAG}_checkm_results_mod2.txt \
    -t ${SLURM_CPUS_PER_TASK} \
    --tab_table

