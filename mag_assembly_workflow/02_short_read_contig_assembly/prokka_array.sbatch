#!/bin/bash

################################################################################
# SCRIPT NAME : prokka_array.sbatch                                            #
# DESCRIPTION : this script annotates and organizes the metagenomic assembly   #
# ARGS        : mapping file (positional, optional)                            #
# AUTHOR      : Matthew C. Hibberd                                             #
# EMAIL       : hibberdm@wustl.edu                                             #
################################################################################

#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -e

if [ -z "${1}" ]; then
    echo "Must pass the metadata"
    exit 1
fi

ml metagenomics_pipeline/3.1

source vars.sh

mkdir -p "$ANNO_DIR"

SID=$( sed -ne '/^#/!p' "${1}" | sed -n "${SLURM_ARRAY_TASK_ID}"p )

ASSEMBLY_FILE="${ASSEMBLY_DIR}/${SID}/${SID}.contigs.fa"
ASSEMBLY_BASE=$( echo "${ASSEMBLY_FILE%.*}" )

#base=${9%.*}
awk '/^>/{print ">" ++i; next}{print}' < "$ASSEMBLY_FILE" > "${ASSEMBLY_BASE}_numbered.fa"

locus_tag="META"$( md5sum "${ASSEMBLY_BASE}_numbered.fa" | awk '{ print toupper($1) }' | cut -c 1-8 )

prefix="NZ_"${SID}

mkdir -p "${ANNO_DIR}/${SID}"

echo -e "${prefix}\t${locus_tag}\n" > "${ANNO_DIR}/${SID}/${prefix}.tag"

if [ -z "${locus_tag}" ]; then
    echo "ERROR: locus_tag not set!"
    exit 1
fi

ml prokka/1.14.0-mod
ml signalp/4.1

prokka \
 --force \
 --outdir "${ANNO_DIR}/${SID}" \
 --prefix "${prefix}" \
 --addgenes \
 --locustag "${locus_tag}" \
 --increment 5 \
 --mincontiglen 500 \
 --genus "Metagenome" \
 --species "metagenome" \
 --strain "${SID}" \
 --cpus "$SLURM_CPUS_PER_TASK" \
 --compliant \
 --centre "WUSM-EFCGSSB" \
 --noanno \
 --notbl2asn \
 "${ASSEMBLY_BASE}_numbered.fa"
