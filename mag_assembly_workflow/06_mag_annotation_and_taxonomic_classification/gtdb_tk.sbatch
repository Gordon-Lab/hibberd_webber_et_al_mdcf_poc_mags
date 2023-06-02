#!/usr/bin/env bash

#SBATCH --mem=200G
#SBATCH --cpus-per-task=24
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-user="hibberdm@wustl.edu"

ml gtdb-tk/1.3.0

export GTDBTK_DATA_PATH="/scratch/ref/gtdbtk/release95/"

gtdbtk \
  classify_wf \
  --genome_dir "/scratch/jglab/hibberdm/metagenomic_assembly/210607_POC_troubleshooting/210609_additional_gtdb_classification/" \
  --out_dir "/scratch/jglab/hibberdm/metagenomic_assembly/210607_POC_troubleshooting/210609_additional_gtdb_classification/" \
  --extension fna \
  --cpus 24
