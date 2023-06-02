
# Gordon Lab Metagenomic Assembly Pipeline
_Author: Matthew C. Hibberd_

_Last update: 2022-03-23_

***Philosophy***: These scripts constitute a workflow for preprocessing, assembly, gene prediction, and quantitation of metagenomic sequencing data. Generally, this workflow follows an "assembly-first" approach, rather than a "mapping-first" approach as is implemented in metaphlan2 and other pipelines. Please note - this workflow is specifically configured for the HTCF computing platform at Washington University in St. Louis School of Medicine and for the specific directory structure used for the analysis described in the manuscript. Additional detail regarding the steps in the pipeline can be found in Extended Data Fig. 1 of Hibberd, Webber, et al.

## Code is archived and maintained at GitLab
[https://gitlab.com/hibberdm/metagenomic_pipeline](https://gitlab.com/hibberdm/metagenomic_pipeline)


## The general steps in the workflow are as follows:

1. **Loading the pipeline and retrieving core workflow components**
	* If you've retrieved this README from the git repo, you can subsequently load the pipeline as below on the HTCF @ WUSTL:
	
	`ml metagenomics_pipeline/<version>`
	
	* alternatively, if you've downloaded the code using git, you'll need to run each step individually using the included sbatch scripts, adapted to your computing infrastructure.
	
	* You can run `retrieve_sbatch.sh` to retrieve example workflow files: 
		* `example_mapping_file.txt` (see next section)
		* `metagenomic_assembly_pipeline.sbatch` (the core workflow script)
		* `config.sh` (a set of workflow script defaults, _see below_)

1. **Mapping file**: Assemble mapping file from sequencing pre/submission information, formatted as in the "example_files" directory in this repo or retrieved via `retrieve_sbatch.sh`. This file contains the following headers (in a commented line) and provides the sample information for all downstream steps in the pipeline. Please **DO NOT** re-order your mapping file columns, as the sbatch scripts that use this information interprets them positionally.

	`#RUN_DIR   ID   SID   SUBGROUP   PLATFORM   HOST`

	Explanation:
	
	* RUN\_DIR => name of sequencing run - more specifically, the folder under raw_data containing the reads for a given run
	* ID => A unique string used to identify pairs of read files corresponding to each sample (usually the index pair string). You do not have to specify multiple lines in the mapping file across the lanes of a single flow cell, but you would for the same sample run across multiple sequencing runs (and/or instruments).
	* SID => an informative SampleID, usually the same used on receipt/documentation of the biological sample.
	* SUBGROUP => [experimental], allows co-assembly of samples in the same subgroup. To assembly singly, use the SID.
	* PLATFORM => either "nextseq", "novaseq", or "generic" - used for file renaming purposes
	* HOST => either "human", "mouse", or "pig" depending on the host your samples were derived from

2. **Configuration**
	* In addition to the mapping file, a number of variables and specifications can be configured in a static file (`vars.sh`) prior to running the pipeline. This file should have been placed in your working directory by running `retrieve_sbatch.sh`.
	* This file contains the set of pipeline defaults for each individual step/script in a series of environmental variables. Any variables/settings specified in this script supersede settings within individual scripts; thus, `vars.sh` serves as a static and portable record of how exactly you ran the workflow and how you deviated from defaults. 
	* The only **REQUIRED** input is `MAP_FILE` (line 17).
	* **OPTIONAL** arguments to be aware of include:
		* `ASSEMBLY_SUBGROUP_FLAG` This boolean enables a branch in the workflow that aggregates data from various SIDs into 'subgroups' as designated by the SUBGROUP column in the mapping file. This choice is useful for MAG assemblies downstream.
		* `ANNO_TBL2ASN_FLAG` This boolean enables use of the NCBI tbl2asn tool in the later stages of the prokka workflow. If you don't need '\<your file>.gff' files or want to save a significant amount of time (for metagenomes) you can set this option to `false` to disable tbl2asn.
		* `ANNO_RRNA_TOOL` This option sets the rRNA prediction tool for prokka. Choices are 'BARRNAP' or 'RNAMMER'. BARRNAP is faster and less license-bound (and purpose-built to integrate with prokka, but RNAMMER is supposedly more accurate, slower, and with more dependencies.

3. **Raw data organization**: 
	* Create a working folder and subfolder to hold the raw data for the analysis 
	
		`mkdir -p <working_dir>/raw_data`
		
	* Copy raw data from `/lts/jglab/` backup into the "raw_data" folder. 

		*Note: If you have multiple runs WITH THE SAME INDEX PAIRS/FILE NAMES, generate subfolders under "raw_data" to hold each run. The pipeline handles the multi-run analysis and concatenation from this point onward.*

4.  **Pipeline initiation**:

* This pipeline is constructed as a series of slurm-managed array jobs, one for each sample in the analysis, to enable maximum parallelization. It uses the "aftercorr" dependency specification in slurm, which enables each sample to proceed through the pipeline steps at a pace in accordance with the complete of each per-sample step.

* The basic pipeline can be run via:

 [**Recommended**] By hardcoding the mapping file into the pipeline script:

		`sbatch metagenomic_pipeline.sbatch`

or, by passing the mapping file into the pipeline script on the command line:

		`sbatch metagenomic_pipeline.sbatch <mapping_file>`

* The pipeline encompasses these steps, each encoded in a separate sbatch script, plus, other separate and optional tools:

	1. Quality trimming and adapter removal: `trim_galore_array.sbatch`
	2. Removal of host-mapping reads: `bowtie_host_filter_array.sbatch`
	3. Concatenation of data across multiple lanes/runs (or subgroups): `cat_multirun_array.sbatch`
	4. Assembly: `megahit_array.sbatch`
	5. Annotation (gene prediction and file organization, not functional predictions by default): `prokka_array.sbatch`
	6. 	Gene quantitation: `bowtie_subread_array.sbatch`
