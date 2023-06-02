#!/bin/bash

##SBATCH --mem=1G

if [ -z "${1}" ]; then
  echo "Must pass the metadata file!"
  exit 1
fi

ml metagenomics_pipeline/3.0

source vars.sh

METADATA=$( sed -ne '/^#/!p' ${1} | sed -n ${SLURM_ARRAY_TASK_ID}p )

IFS=$'\t' read -r RUN_DIR ID SID SUBGROUP PLATFORM HOST <<< "$METADATA"
unset IFS

echo -e "ID\t$ID"
echo -e "Platform\t$PLATFORM"

# Example filenames
#Gordon_9B_N701_S613_TAAGGCGA_GGAGTTCC_S33_L001_R1_001.fastq.gz
#Gordon_Carrie001_N709_S508_GCTACGCT_AGGCTTAG_S2_R1_001.fastq.gz
#Gordon_Hannah045_N821_S536_CACGAGTT_TTAACCGA_S46_R1_001.fastq

if [[ "$PLATFORM" == "nextseq" ]]; then

  echo "NextSeq run : Using CGS naming convention."

  R1_count=$( ls ${READS_DIR}/${RUN_DIR}/*${ID}*R1_001.fastq* | wc -l )
  R2_count=$( ls ${READS_DIR}/${RUN_DIR}/*${ID}*R2_001.fastq* | wc -l )

  #    file s_4_4_withindex_sequence.txt_TTCTTCTA-TGTAATGA.fq.gz | grep -c  "gzip"

  if [[ "$R1_count" = 1 && "$R2_count" = 1 ]]; then


    R1=$( find ${READS_DIR}/${RUN_DIR}/ -type f -name "*${ID}*R1_001.fastq*" )
    R2=$( find ${READS_DIR}/${RUN_DIR}/ -type f -name "*${ID}*R2_001.fastq*" )


    if [ -f ${READS_DIR}/${RUN_DIR}/s_1_1_withindex_sequence.txt_${ID}.fq.gz ]; then
      echo "${READS_DIR}/${RUN_DIR}/s_1_1_withindex_sequence.txt_${ID}.fq.gz  exists, no need to rename"
    else
      echo "Linking ${R1} to ${READS_DIR}/${RUN_DIR}/s_1_1_withindex_sequence.txt_${ID}.fq.gz"
      ln -s ${R1} ${READS_DIR}/${RUN_DIR}/s_1_1_withindex_sequence.txt_${ID}.fq.gz
    fi

    if [ -f ${READS_DIR}/${RUN_DIR}/s_1_4_withindex_sequence.txt_${ID}.fq.gz ]; then
      echo "${READS_DIR}/${RUN_DIR}/s_1_4_withindex_sequence.txt_${ID}.fq.gz exissts, no need to rename"
    else
      echo "Linking ${R2} to ${READS_DIR}/${RUN_DIR}/s_1_4_withindex_sequence.txt_${ID}.fq.gz"
      ln -s ${R2} ${READS_DIR}/${RUN_DIR}/s_1_4_withindex_sequence.txt_${ID}.fq.gz
    fi

  else

    if [[ "$R1_count" > 1 || "$R2_count" > 1 ]]; then

      echo "[ERROR] More than one reads file matching ${ID} R1 or R2 found! Please confirm pattern match"
      exit 1

    else

      echo "[ERROR] No reads file matching ${ID} R1 or R2 found! Please confirm pattern match"
      exit 1

    fi

  fi

elif [[ "$PLATFORM" == "novaseq" ]]; then

  echo "NovaSeq run : New filename format"
  # Example filename : AAGAGGCA-ACTGCATA_S2_L003_R1_001.fastq.gz or 78149-9_TAAGGCGA-CTCTCTAT_S664_L001_R2_001.fastq.gz"
  ## Loop through all lanes

  for ((i = 0; i < 4; ++i));
  do
    LANE=$(( i + 1 ))

    if [[ -f $( find ${READS_DIR}/${RUN_DIR}/ -type f -name "*${ID}*L00${LANE}*R1_001.fastq*") &&
    -f $( find ${READS_DIR}/${RUN_DIR}/ -type f -name "*${ID}*L00${LANE}*R2_001.fastq*") ]]; then
      #if [[ -f "${READS_DIR}/${RUN_DIR}/${ID}_S??_L00${LANE}_R1_001.fastq.gz" &&
      #-f "${READS_DIR}/${RUN_DIR}/${ID}*"L00"${LANE}"_R2_001.fastq.gz" ]]; then

      R1_count=$( ls ${READS_DIR}/${RUN_DIR}/*${ID}*L00${LANE}*R1_001.fastq* | wc -l )
      R2_count=$( ls ${READS_DIR}/${RUN_DIR}/*${ID}*L00${LANE}*R2_001.fastq* | wc -l )

      if [[ "$R1_count" = 1 && "$R2_count" = 1 ]]; then

        R1=$( find ${READS_DIR}/${RUN_DIR}/ -type f -name "*${ID}*L00${LANE}*R1_001.fastq*" )
        R2=$( find ${READS_DIR}/${RUN_DIR}/ -type f -name "*${ID}*L00${LANE}*R2_001.fastq*" )


        if [ -f ${READS_DIR}/${RUN_DIR}/s_${LANE}_1_withindex_sequence.txt_${ID}.fq.gz ]; then
          echo "${READS_DIR}/${RUN_DIR}/s_${LANE}_1_withindex_sequence.txt_${ID}.fq.gz  exists, no need to rename"
        else
          echo "Linking ${R1} to ${READS_DIR}/${RUN_DIR}/s_${LANE}_1_withindex_sequence.txt_${ID}.fq.gz"
          ln -s ${R1} ${READS_DIR}/${RUN_DIR}/s_${LANE}_1_withindex_sequence.txt_${ID}.fq.gz
          #echo "${R1} ${READS_DIR}/${RUN_DIR}/s_${LANE}_1_withindex_sequence.txt_${ID}.fq.gz"
        fi

        if [ -f ${READS_DIR}/${RUN_DIR}/s_${LANE}_4_withindex_sequence.txt_${ID}.fq.gz ]; then
          echo "${READS_DIR}/${RUN_DIR}/s_${LANE}_4_withindex_sequence.txt_${ID}.fq.gz exists, no need to rename"
        else
          echo "Linking ${R2} to ${READS_DIR}/${RUN_DIR}/s_${LANE}_4_withindex_sequence.txt_${ID}.fq.gz"
          ln -s ${R2} ${READS_DIR}/${RUN_DIR}/s_${LANE}_4_withindex_sequence.txt_${ID}.fq.gz
          #echo "ln -s ${R2} ${READS_DIR}/${RUN_DIR}/s_${LANE}_4_withindex_sequence.txt_${ID}.fq.gz"
        fi

      else

        if [[ "$R1_count" > 1 || "$R2_count" > 1 ]]; then
          echo "[ERROR] More than one reads file matching ${ID} R1 or R2 found! Please confirm pattern match"
          exit 1
        else
          echo "[ERROR] No reads file matching ${ID} R1 or R2 found! Please confirm pattern match"
          exit 1

        fi
      fi
    fi
  done

elif [[ "$PLATFORM" == "generic" ]]; then

  echo "Generic run : Using ID_R1.fastq formatting"
  # Example filename : AAGAGGCA-ACTGCATA_S2_L003_R1_001.fastq.gz

  ## Loop through all lanes

  if [[ -f $( find ${READS_DIR}/${RUN_DIR}/ -type f -name "${ID}_R1.fastq*") &&
  -f $( find ${READS_DIR}/${RUN_DIR}/ -type f -name "${ID}_R2.fastq*") ]]; then

    R1_count=$( ls ${READS_DIR}/${RUN_DIR}/${ID}_R1.fastq* | wc -l )
    R2_count=$( ls ${READS_DIR}/${RUN_DIR}/${ID}_R2.fastq* | wc -l )

    if [[ "$R1_count" = 1 && "$R2_count" = 1 ]]; then

      R1=$( find ${READS_DIR}/${RUN_DIR}/ -type f -name "${ID}_R1.fastq*" )
      R2=$( find ${READS_DIR}/${RUN_DIR}/ -type f -name "${ID}_R2.fastq*" )

      if [ -f ${READS_DIR}/${RUN_DIR}/s_1_1_withindex_sequence.txt_${ID}.fq.gz ]; then
        echo "${READS_DIR}/${RUN_DIR}/s_1_1_withindex_sequence.txt_${ID}.fq.gz  exists, no need to rename"
      else
        echo "Linking ${R1} to ${READS_DIR}/${RUN_DIR}/s_1_1_withindex_sequence.txt_${ID}.fq.gz"
        ln -s ${R1} ${READS_DIR}/${RUN_DIR}/s_1_1_withindex_sequence.txt_${ID}.fq.gz
        #echo "${R1} ${READS_DIR}/${RUN_DIR}/s_${LANE}_1_withindex_sequence.txt_${ID}.fq.gz"
      fi

      if [ -f ${READS_DIR}/${RUN_DIR}/s_1_4_withindex_sequence.txt_${ID}.fq.gz ]; then
        echo "${READS_DIR}/${RUN_DIR}/s_1_4_withindex_sequence.txt_${ID}.fq.gz exists, no need to rename"
      else
        echo "Linking ${R2} to ${READS_DIR}/${RUN_DIR}/s_1_4_withindex_sequence.txt_${ID}.fq.gz"
        ln -s ${R2} ${READS_DIR}/${RUN_DIR}/s_1_4_withindex_sequence.txt_${ID}.fq.gz
        #echo "ln -s ${R2} ${READS_DIR}/${RUN_DIR}/s_${LANE}_4_withindex_sequence.txt_${ID}.fq.gz"
      fi

    else

      if [[ "$R1_count" > 1 || "$R2_count" > 1 ]]; then
        echo "[ERROR] More than one reads file matching ${ID} R1 or R2 found! Please confirm pattern match"
        exit 1
      else
        echo "[ERROR] No reads file matching ${ID} R1 or R2 found! Please confirm pattern match"
        exit 1

      fi
    fi
  fi  
fi
