#!/usr/bin/env bash

################################################################################
# SCRIPT NAME : kraken_tax.sbatch                                              #
# DESCRIPTION : this script performs a taxonomic analysis of sequencing reads  #
# ARGS        : mapping file (positional)                                      #
# AUTHOR      : Matthew C. Hibberd                                             #
# EMAIL       : hibberdm@wustl.edu                                             #
################################################################################

#SBATCH --mem=1G

set -e

if [ -z "${1}" ]; then
    echo "Must pass the metadata"
    exit 1
fi

ml metagenomics_pipeline/3.0
ml prodigal

source vars.sh

MAG=$( sed -ne '/^#/!p' "${1}" | sed -n "${SLURM_ARRAY_TASK_ID}"p )

echo -e "MAG: ${MAG}"

PROD_DIR=${ANNO_DIR}/prodigal_partial/

mkdir -p $PROD_DIR

prodigal \
  -i ${ANNO_DIR}/${MAG}/${MAG}.fna \
  -m \
  -g 11 \
  -p single \
  -f gff \
  -q \
  > ${PROD_DIR}/${MAG}_prodigal.gff

prodigal \
  -i ${ANNO_DIR}/${MAG}/${MAG}.fna \
  -m \
  -g 11 \
  -p single \
  -f gbk \
  -q \
  > ${PROD_DIR}/${MAG}_prodigal.gbk
