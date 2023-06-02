# Code relevant to Hibberd, Webber et al: "_Bioactive glycans in a microbiome-directed food for malnourished children_"
_Author: Matthew C. Hibberd_

_Last update: 2023-06-02_

***Preamble***: We present the include R, python and bash code as guidelines and examples of the bioinformatics analyses we performed to generate MAGs, quantify their DNA and transcript abundances, plus analyze these datasets in the context of features of their genomes and the anthropometric characteristics of their hosts. We also include an R codebase detailing analyses related to glycan abundances in diet and feces. With the included bash code, sequencing datasets obtained from ENA (PRJEB45356) and reasonable knowledge of high-performance computing, the software included can be adapted to replicate our sequencing-based DNA and RNA analyses. With the included R code and abundance and metadata datasets from our Supplementary Materials, one should be able to replicate the various statistical analyses and plots generated for the manuscript. Most individual statistical analyses can be run in <1 hr on typical desktop/laptop hardware. Steps related to MAG assembly, metatranscriptomic analyses, etc. can be expected to take ~1-24 hours, approximately, depending on parallelization ability and available processors.

## Code is archived and maintained at GitLab:
[https://gitlab.com/hibberdm/hibberd\_webber\_et\_al\_mdcf\_poc\_mags](https://gitlab.com/hibberdm/hibberd_webber_et_al_mdcf_poc_mags)


## The included code (in folders) and their use is as follows:

1. **anthropometry\_analysis**
	* R code relevant for calculating beta-WLZ for each child, plus apportioning children into response quartiles based on their WLZ response.
2. **cami\_mag\_quantitation\_benchmarking**
	* R code, Python code, and bash code related to defining and supporting the contig and MAG quantitation strategy defined in Extended Data Figs. 1 and 11, plus Supplementary Discussion.
3. **diet\_intake**
	* Data and R code related to analysis of the food frequency questionnaire.
4. **mag\_abundance\_analysis**
	* R code related to analysis of MAG abundances versus anthropometry, plus predicted phenotype enrichments and data wrangling for generating supplementary tables
5. **mag\_assembly\_workflow**
	* Bash and python code related to the Gordon Lab MAG assembly workflow. This folder has its own README.
6. **metatranscriptomics\_workflow**
	* Bash and R code related to transcript quantitation from RNA-sequencing data, plus code relevant to analysis of differential expression between diet treatment groups and response quartiles.
7. **monosaccharide\_and\_glycosidic\_linkage\_analysis**
	* Data and R code related to analysis of MS-based linkage, mono- and polysaccharide quantitation in diets and in the feces of trial participants.
8. **prevotella\_copri\_puls**
	* R code related to plotting PUL conservation profiles for P. copri MAGs and isolates, plus code relevant for relating transcriptional levels of PUL transcripts to other data types.

## Software versions and installation:
* Software versions for R and various HTCF-relevant software packages are stated in the Methods of Hibberd, Webber et al. and in the included scripts. Installation will depend on each user's specific hardware, but is easily accomplished by standard routes. Custom R scripts can be run without installation after installing the included libraries/packages as detailed in each script and Methods. R installation details for various operating systems can be found at [https://cran.r-project.org/](https://cran.r-project.org/). Details of the WUSM High-Throughput Computing Facility configuration can be found at [https://htcf.wustl.edu/docs/](https://htcf.wustl.edu/docs/).