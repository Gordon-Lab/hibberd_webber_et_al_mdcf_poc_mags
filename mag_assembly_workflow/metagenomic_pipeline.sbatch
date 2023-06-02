#!/bin/bash

################################################################################
# SCRIPT NAME : metagenomic_pipeline.sbatch                                    #
# DESCRIPTION : pipeline submission template for metagenomic sequence analysis #
# ARGS        : mapping file (positional, optional)                            #
# AUTHOR      : Matthew C. Hibberd                                             #
# EMAIL       : hibberdm@wustl.edu                                             #
################################################################################

#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-user="$USER@wustl.edu"

## SHELL BEHAVIOR MODIFIERS
# 'set -e' causes immediate exit on non-zero status
set -e

## LOAD THE MODULE
ml metagenomics_pipeline/3.1

## VARIABLE HANDLING
source vars.sh

## DIRECTORY INITIALIZATION
# Log directory is created up front, other directories are created as needed
mkdir -p $LOG_DIR

## USER INPUT
# Enter mapping file name (if in WORK_DIR) or full path OR read in from command line
MAP_FILE=""
# Choose "SAMPLE" or "SUBGROUP"
ASSEMBLY_STRATEGY="SAMPLE"

## MAPPING FILE QC
# Make sure the mapping file exists
if [ -z "${MAP_FILE}" ]; then
    MAP_FILE=$1
    if [ -z "${MAP_FILE}" ]; then
        echo "ERROR: Must specify or pass the mapping file!"
        exit 1
    fi
fi

# Convert newlines from dos to unix format (if necessary)
dos2unix -ic ${MAP_FILE} | xargs dos2unix

# This command makes sure a newline exists at the end of the mapping file.
sed -i -e '$a\' ${MAP_FILE}

## JOB ARRAY HANDLING
# This command determines the number of samples from the mapping file for use in
# creating job arrays.
SAMPLES=$( sed -ne '/^#/!p' "${MAP_FILE}" | wc -l )

if [ -f "${WORK_DIR}/SID_list.txt" ]; then
    
    if [ -n "$(grep "^[^#]" ${MAP_FILE} | cut -f 3 | sort | uniq | cmp "${WORK_DIR}/SID_list.txt" - )" ]; then
      echo "[STATUS] Existing SID_list.txt doesn't match. Archiving existing and creating a new SID_list.txt: ${WORK_DIR}/SID_list.txt"
      mv "${WORK_DIR}/SID_list.txt" "${WORK_DIR}/archive_SID_list.txt"
      grep "^[^#]" ${MAP_FILE} | cut -f 3 | sort | uniq > "${WORK_DIR}/SID_list.txt"
    else
      echo "[STATUS] Using existing SID_list.txt: ${WORK_DIR}/SID_list.txt"
    fi
else
    echo "[STATUS] Creating new SID_list.txt: ${WORK_DIR}/SID_list.txt"
    grep "^[^#]" ${MAP_FILE} | cut -f 3 | sort | uniq > "${WORK_DIR}/SID_list.txt"
fi

ASSEMBLIES=$( wc -l "${WORK_DIR}/SID_list.txt" | cut -d " " -f 1 )

## PREPROCESSING
J1=$( sbatch -o "${LOG_DIR}/slurm-nextseq_rename-%A_%a.out" --array=1-"${SAMPLES}"%20 nextseq_renamer.sbatch "${MAP_FILE}" | awk '{ print $NF }' )

J2=$( sbatch -o "${LOG_DIR}/slurm-trim_galore-%A_%a.out" --array=1-"${SAMPLES}"%100 -d aftercorr:"${J1}" trim_galore_array.sbatch "${MAP_FILE}" | awk '{ print $NF }' )
J3=$( sbatch -o "${LOG_DIR}/slurm-bowtie_host_filter-%A_%a.out" --array=1-"${SAMPLES}"%100 -d aftercorr:"${J2}" bowtie_host_filter_array.sbatch "${MAP_FILE}" | awk '{ print $NF }' )
J4=$( sbatch -o "${LOG_DIR}/slurm-cat-%A_%a.out" --array=1-"${ASSEMBLIES}"%100 -d afterok:"${J3}" cat_multirun_array.sbatch "${WORK_DIR}/SID_list.txt" | awk '{ print $NF }' )

## ASSEMBLY, ANNOTATION, AND COUNTING
J5=$( sbatch -o "${LOG_DIR}/slurm-megahit-%A_%a.out" --array=1-"${ASSEMBLIES}"%100 -d aftercorr:"${J4}" megahit_array.sbatch "${WORK_DIR}/SID_list.txt" | awk '{ print $NF }' )
J6=$( sbatch -o "${LOG_DIR}/slurm-prokka-%A_%a.out" --array=1-"${ASSEMBLIES}"%100 -d aftercorr:"${J5}" prokka_array.sbatch "${WORK_DIR}/SID_list.txt" | awk '{ print $NF }' )
J7=$( sbatch -o "${LOG_DIR}/slurm-subread-%A_%a.out" --array=1-"${ASSEMBLIES}"%100 -d aftercorr:"${J6}" bowtie_subread_array.sbatch "${WORK_DIR}/SID_list.txt" | awk '{ print $NF }' )

## OPTIONAL STEPS
#sbatch -o "${LOG_DIR}/slurm-tax_motus2-%A_%a.out" --array=1-"${ASSEMBLIES}"%100 -d aftercorr:"${J5}" motus_v2.sbatch "${WORK_DIR}/SID_list.txt"
#sbatch -o "${LOG_DIR}/slurm-tax_metaphlan2-%A_%a.out" --array=1-"${ASSEMBLIES}"%100 -d aftercorr:"${J5}" metaphlan2.sbatch "${WORK_DIR}/SID_list.txt"
#sbatch -o "${LOG_DIR}/slurm-tax_kraken2-%A_%a.out" --array=1-"${ASSEMBLIES}"%100 -d aftercorr:"${J5}" kraken2_tax.sbatch "${WORK_DIR}/SID_list.txt"
#sbatch -o "${LOG_DIR}/slurm-plass-%A_%a.out" --array=1-"${ASSEMBLIES}"%100 -d aftercorr:"${J5}" plass.sbatch "${WORK_DIR}/SID_list.txt"

## QC
#sbatch -o "${LOG_DIR}/slurm-nonpareil-%A_%a.out" --array=1-"${ASSEMBLIES}"%100 -d aftercorr:"${J4}" nonpareil.sbatch "${WORK_DIR}/SID_list.txt"
#sbatch -o logs/slurm-metaquast.out ~/metagenomic_pipeline/bin/metaquast.sbatch

#mcSEED
#KEGG
#QC - read counting, % mapping to assemblies, etc.
#Metaquast
#Cleanup

## MISCELLANEOUS
#-size +0
#find . -size +0 -type f -name 'contig.fa' | sed -r 's|/[^/]+$||' |sort |uniq | wc -l
