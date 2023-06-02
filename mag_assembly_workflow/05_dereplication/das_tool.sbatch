#!/bin/bash

#Slurm email messages
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-user="hibberdm@wustl.edu"
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

ml metagenomics_pipeline/3.1
ml R/3.6.1
ml pullseq/1.0.2
ml prodigal/2.6.3
#ml usearch/11.0.667
ml diamond/0.9.34-python-3.6.5

source vars.sh

BIN_DIR_METABAT="09_contig_binning_metabat"
BIN_DIR_MAXBIN="09_contig_binning_maxbin"

DAS_HOME="/scratch/jglab/hibberdm/DAS_Tool"
DAS_OUT="10_DAS"

SUBGROUP=$( sed -ne '/^#/!p' "${1}" | sed -n "${SLURM_ARRAY_TASK_ID}"p )

mkdir -p $DAS_OUT
mkdir -p ${DAS_OUT}/NZ_${SUBGROUP}

${DAS_HOME}/src/Fasta_to_Scaffolds2Bin.sh \
  -e fa \
  -i ${BIN_DIR_METABAT}/NZ_${SUBGROUP}/ \
  > ${DAS_OUT}/NZ_${SUBGROUP}/NZ_${SUBGROUP}_metabat.s2b

${DAS_HOME}/src/Fasta_to_Scaffolds2Bin.sh \
  -e fasta \
  -i ${BIN_DIR_MAXBIN}/NZ_${SUBGROUP}/ \
  > ${DAS_OUT}/NZ_${SUBGROUP}/NZ_${SUBGROUP}_maxbin.s2b

${DAS_HOME}/DAS_Tool \
  -i ${DAS_OUT}/NZ_${SUBGROUP}/NZ_${SUBGROUP}_metabat.s2b,${DAS_OUT}/NZ_${SUBGROUP}/NZ_${SUBGROUP}_maxbin.s2b \
  -c ${ANNO_DIR}/${SUBGROUP}/NZ_${SUBGROUP}.fna \
  -l NZ_${SUBGROUP}_metabat,NZ_${SUBGROUP}_maxbin \
  --write_bin_evals 1 \
  --write_bins 1 \
  --db_directory $DAS_HOME \
  -t ${SLURM_CPUS_PER_TASK} \
  -o ${DAS_OUT}/NZ_${SUBGROUP}/NZ_${SUBGROUP} \
  --search_engine diamond

