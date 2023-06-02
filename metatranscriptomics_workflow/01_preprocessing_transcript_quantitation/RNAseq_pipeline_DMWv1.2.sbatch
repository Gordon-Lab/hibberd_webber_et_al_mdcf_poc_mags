#!/usr/bin/env bash

################################################################################
# SCRIPT NAME : RNAseq_pipeline_DMWv1.2.sbatch                                 #
# DESCRIPTION : RNAseq pipeline w/ kallisto pseudoalignment                    #
# AUTHOR      : Daniel M Webber                                                #
# EMAIL       : dmwebber@wustl.edu                                             #
################################################################################

#SBATCH --array=1-116%116   #Create as many tasks as lines in the RawSequences_list (Example: --array=1-116%116)
#SBATCH --cpus-per-task=5
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-user="user@wustl.edu"
#SBATCH --mem=100G
####SBATCH --mem-per-cpu=20G
#SBATCH -o /path/to/directory/logfile_%A_%a.txt # STDOUT

###################
## Setup          #
###################
# Make the kallisto index before running this cript. General instructions for making the kallisto index are below
#   srun --mem=50G --cpus-per-task=16 -J interactive -p interactive --pty /bin/bash
#   concatenated_MAGs="/path/to/directory/MAGs_build1.ffn"	
#   kallisto_index_name="/path/to/directory/index.kallisto"
#   ml kallisto/0.43.0   
#   kallisto index -i ${kallisto_index_name} ${concatenated_MAGs}
	
###################
## 1. User Input ##
###################
sample_list="/path/to/directory/SampleSheet.csv"	# Supply a csv list of input files with R1,R2 in columns 1 and 2 (use absolute paths, files can be either gz compressed or uncompressed)
kallisto_index_name="/path/to/directory/index.kallisto"	# Name for the kallisto index
kallisto_out_dir="/path/to/directory/index.kallisto"     # This directory will be made if it does not already exist
rRNA_tags="/path/to/directory/index.kallisto/ribosomal_RNA_allMAGs.txt"    # Supply a line-delimited list of locus tags that designate rRNA  
QC_out_dir="/path/to/directory/index.kallisto/QC/"

######################
## 2. File handling ##
######################     
## 2a. Parse the file list and make directories as needed
lineN=$( sed -n ${SLURM_ARRAY_TASK_ID}p $sample_list )
R1=`echo $lineN | awk -F "," '{ print $1 }'`
R2=`echo $lineN | awk -F "," '{ print $2 }'`
R1_filename=`basename ${R1}`
R2_filename=`basename ${R2}`
SID=`echo $R1_filename | awk -F "_R1" '{ print $1 }'`
mkdir ${kallisto_out_dir}
mkdir ${QC_out_dir}
QC_out_dir_folder=${QC_out_dir}/${SID}
mkdir ${QC_out_dir_folder}
kallisto_out_dir_folder=${kallisto_out_dir}/${SID}
mkdir ${kallisto_out_dir_folder}
QCfolder=${QC_out_dir_folder}/QC
mkdir ${QCfolder}


## 2b. Print directory info to std out
Nlines=`wc -l < ${sample_list}`
echo 'Lines in sample file:' ${Nlines}
echo 'Working on task' ${SLURM_ARRAY_TASK_ID} #Update to STDOUT
echo 'Input directory:' ${dir}
echo 'SID:' ${SID}
echo 'R1:' ${R1}
echo 'R2:' ${R2}
echo 'Kallisto index name:' ${kallisto_index_name} 
echo 'Sample output folder:' ${QC_out_dir_folder}
echo 'Sample QC folder:' ${QCfolder}

#############
### 3.  QC ##
#############
### 3a. Run FastQC on raw sequence data
echo 'Run FastQC on raw-sequence data'
trim_report=${QCfolder}/${SID}_R1.fastq.gz_trimming_report.txt

if [ -f "$R1trim" ]; then 
	echo 'Previously performed, skipping this step'  
else
	ml fastqc/0.11.7-java-11
	fastqc ${R1} --outdir=${QCfolder}
	fastqc ${R2} --outdir=${QCfolder}
fi

### 3b. Trim indexes and adapers with trim galore
echo 'Trim indexes and adapters with Trim galore'

#Make variables to hold the full PathName for trimmed R1 and R2 fastq sequences
R1trim=${QCfolder}/${SID}_R1_val_1.fq.gz
R2trim=${QCfolder}/${SID}_R2_val_2.fq.gz

if [ -f "$R1trim" ]; then
	echo 'Previously performed, skipping this step'
else 
	ml trim_galore/0.6.4
	trim_galore \
	--cores "${SLURM_CPUS_PER_TASK}" \
	--phred33 \
	--length 100 \
	--2colour 20 \
	--paired \
	--basename "${SID}" \
	--output_dir "${QCfolder}" \
	--trim-n \
	"${R1}" \
	"${R2}"
fi

### 3c. Run FastQC on sequence data
echo 'Run FastQC on trimmed-sequence data'

if [ -f "${QCfolder}/${SID}_R1_val_1_fastqc.html" ]; then
	echo 'Previously performed, skipping this step'
else
	ml fastqc/0.11.7-java-11
	fastqc ${R1trim} --outdir=${QCfolder}
	fastqc ${R2trim} --outdir=${QCfolder}
fi

#####################################
### 4. Pseudocounts with kallisto  ##
#####################################
### 4a. Quantification with kallisto
ml kallisto/0.43.0

echo 'Determine pseudocount using kallisto'

kallisto quant \
-i ${kallisto_index_name} \
-o ${kallisto_out_dir_folder} \
--bootstrap-samples=1 \
--rf-stranded \
--threads=${SLURM_CPUS_PER_TASK} \
${R1trim} \
${R2trim}

### 4b. Use grep to filter out counts from ribosomal_RNA 
#Add headers
printf "target_id\tlength\teff_length\test_counts\ttpm\n" > ${kallisto_out_dir_folder}/abundance_rRNA_filtered.tsv
printf "target_id\tlength\teff_length\test_counts\ttpm\n" > ${kallisto_out_dir_folder}/abundance_rRNA.tsv

#Grep
grep -vf ${rRNA_tags} ${kallisto_out_dir_folder}/abundance.tsv >> ${kallisto_out_dir_folder}/abundance_rRNA_filtered.tsv 
grep -f ${rRNA_tags} ${kallisto_out_dir_folder}/abundance.tsv >> ${kallisto_out_dir_folder}/abundance_rRNA.tsv

### 4c. Add up total mapped reads, total reads not mapped to tRNA, tRNA mapped reads:
non_rRNA_count=`cat ${kallisto_out_dir_folder}/abundance_rRNA_filtered.tsv | awk 'NR > 1 { sum += $4 } END { printf "%d\n", sum }'`
rRNA_count=`cat ${kallisto_out_dir_folder}/abundance_rRNA.tsv | awk 'NR > 1 { sum += $4 } END { printf "%d\n", sum }'`
total_read_count=`cat ${kallisto_out_dir_folder}/abundance.tsv | awk 'NR > 1 { sum += $4 } END { printf "%d\n", sum }'`

####################################
# Output the read count:
####################################
pre_trim_ct=`bc <<< $(zcat $R1|wc -l)/4`
post_trim_ct=`bc <<< $(zcat $R1trim|wc -l)/4`

# append this info to a single file causes problems with parrell processing so make individual files and cat
printf "$SID\t$pre_trim_ct\t$post_trim_ct\t$total_read_count\t$non_rRNA_count\t$rRNA_count\n" > "${kallisto_out_dir}/${SID}_temp_counts.tsv"

### Once reaching the last sample, then compile the read count info:
if [[ ${SLURM_ARRAY_TASK_ID} -eq ${Nlines} ]]
then
	echo 'Wait to make sure that all processes are complete'
	sleep 40m
	echo 'concatenating counts'	
	printf "Sample\tCount_pre_filter\tCount_post_filter\tTotal_mapped_reads\tNon_rRNA_count\trRNA_count\n" > ${kallisto_out_dir}/read_ct.tsv
	cat ${kallisto_out_dir}/*_temp_counts.tsv >> ${kallisto_out_dir}/read_ct.tsv
	sleep 1m
	#rm ${kallisto_out_dir}/*_temp_counts.tsv
	echo 'Done with QC and kallisto read counts for this sample'
 else
	echo 'Done with QC and kallisto read counts for this sample'
 fi
