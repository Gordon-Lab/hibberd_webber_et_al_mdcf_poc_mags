#!/usr/bin/env bash

################################################################################
# SCRIPT NAME : kraken_tax.sbatch                                              #
# DESCRIPTION : this script performs a taxonomic analysis of sequencing reads  #
# ARGS        : mapping file (positional)                                      #
# AUTHOR      : Matthew C. Hibberd                                             #
# EMAIL       : hibberdm@wustl.edu                                             #
################################################################################

#SBATCH --mem=200G
#SBATCH --cpus-per-task=12

set -e

if [ -z "${1}" ]; then
    echo "Must pass the metadata"
    exit 1
fi

ml metagenomics_pipeline/3.0
ml kraken2/2.0.8-beta
ml bracken/2.5

source vars.sh

mkdir -p "$TAX_DIR"

MAG=$( sed -ne '/^#/!p' "${1}" | sed -n "${SLURM_ARRAY_TASK_ID}"p )

echo -e "MAG: ${MAG}"

KRAKEN_DB="/scratch/jglab/hibberdm/kraken2_db/kraken/"
BRACKEN_DB="/scratch/jglab/hibberdm/kraken2_db/bracken/"
ANNO_DIR="13_bin_annotation_1000"

#kraken2-inspect --db $KRAKEN_DB

export TAX_DIR="1000_bin_kraken2_taxonomy"

mkdir -p $TAX_DIR

kraken2 \
    --db $KRAKEN_DB \
    --threads $SLURM_CPUS_PER_TASK \
    --report ${TAX_DIR}/${MAG}.k2report \
    --use-names \
    ${ANNO_DIR}/${MAG}/${MAG}.fna \
    --confidence 0.1 \
    --output ${TAX_DIR}/${MAG}.kraken2

bracken \
    -d $BRACKEN_DB \
    -i ${TAX_DIR}/${MAG}.k2report \
    -o ${TAX_DIR}/${MAG}.bracken

