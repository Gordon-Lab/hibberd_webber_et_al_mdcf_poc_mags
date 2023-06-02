#!/bin/bash
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-user="hmbucklin@wustl.edu"

mkdir -p prokka_annotations

SID="CAMI_data"

ASSEMBLY_FILE="/scratch/jglab/hmbucklin/03032021_CAMI_mouse_reanalysis/791_genomes_w_Ns/all_791_genomes_w_Ns_cat.fna"
ASSEMBLY_BASE=$( echo "${ASSEMBLY_FILE%.*}" )

#base=${9%.*}
awk '/^>/{print ">" ++i; next}{print}' < "$ASSEMBLY_FILE" > "${ASSEMBLY_BASE}_numbered.fa"

locus_tag="META"$( md5sum "${ASSEMBLY_BASE}_numbered.fa" | awk '{ print toupper($1) }' | cut -c 1-8 )

prefix="CAMI"

if [ -z "${locus_tag}" ]; then
    echo "ERROR: locus_tag not set!"
    exit 1
fi

ml prokka/1.14.0-mod
ml signalp/4.1

prokka \
 --force \
 --outdir "/scratch/jglab/hmbucklin/03032021_CAMI_mouse_reanalysis/prokka_annotations" \
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
 /scratch/jglab/hmbucklin/03032021_CAMI_mouse_reanalysis/791_genomes_w_Ns/all_791_genomes_w_Ns_cat.fna
 