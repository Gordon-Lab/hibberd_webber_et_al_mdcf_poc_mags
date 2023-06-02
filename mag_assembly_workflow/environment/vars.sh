#!/bin/sh

## DIRECTORIES

WORK_DIR=$( pwd )
export WORK_DIR
export LOG_DIR=${WORK_DIR}/logs
export READS_DIR=${WORK_DIR}/raw_data
export QC_DIR=${WORK_DIR}/01_quality_control
export HF_DIR=${WORK_DIR}/02_host_filtered
export ASSEMBLY_DIR=${WORK_DIR}/03_assembly
export ANNO_DIR=${WORK_DIR}/04_annotation
export ADV_ANNO_DIR=${WORK_DIR}/05_advanced_annotation
export COUNTS_DIR=${WORK_DIR}/06_counts
export TAX_DIR=${WORK_DIR}/07_taxonomic_profiling
export CLUSTER_DIR=${WORK_DIR}/08_protein_clustering
export BIN_DIR=${WORK_DIR}/09_contig_binning

## OTHER VARIABLES

export HOST_DB_HUMAN="/scratch/ref/genomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
export HOST_DB_MOUSE="/scratch/ref/genomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome"
export HOST_DB_PIG="/scratch/ref/genomes/Sus_scrofa/UCSC/susScr3/Sequence/Bowtie2Index/genome"

echo "[VARIABLES] Environment variables set."
