library(fgsea)
library(tidyverse)
library(readxl)
library(progress)
library(ggrepel)
#library(ggtree)
#library(treeio)
library(progress)
library(DESeq2)
library(edgeR)

setwd("~/Documents/Projects/human_studies/MDCF POC/210608_bin_abundance_agg_and_calc_build4/")

# set useful variables for file naming and version tracking
datestring <- format(Sys.time(), "%y%m%d_%H%M")
study <- "MDCF_POC"
experiment <- "wlz-MAGs"

# custom, colorblind-friendly plot color palette
cbbPalette <- c("#999999", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "darkred", "darkblue", "darkgrey", "white")

# pre-built functions stored in an external file for extracting information from linear model objects
source("functions.R")

# reading in metadata and abundance data ----------------------------------

# read in mapping file/metadata for the samples
map <- read_excel("~/Documents/Projects/human_studies/MDCF POC/200206 16S ASV Analysis/200428_MDCF_POC_fecal_plasma_metadata.xlsx") %>%
  filter(sample_type == "fecal")

# read in anthropometry data
map.anthropometry <- read.table("~/Documents/Projects/human_studies/MDCF POC/200423_MDCF_POC_anthropometry.txt", sep = "\t", header = T)

# read in taxonomic classifications for all MAGs from GTDB-tk. These were performed in two batches, so we'll combine the files
tax1 <- read.table("~/Documents/Projects/human_studies/MDCF POC/210115_kallisto_mag_purity/gtdb_checkm_magpurify_clean/gtdbtk.bac120.summary.tsv", sep = "\t", header = TRUE) %>%
  separate(classification, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
  mutate(bin = gsub("^NZ_", "", user_genome))
tax2 <- read.table("~/Documents/Projects/human_studies/MDCF POC/210608_bin_abundance_agg_and_calc_build4/210610_gtdb_tk_33_assignments/gtdbtk.bac120.summary.tsv", sep = "\t", header = TRUE) %>%
  separate(classification, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
  mutate(bin = gsub("^NZ_", "", user_genome))
tax <- bind_rows(tax1, tax2)

# identify the individual, per-sample counts files from kallisto.
files <- list.files(path = "kallisto_counts_1000", 
                    pattern = "abundance.tsv", 
                    full.names = TRUE, 
                    recursive = TRUE)

# each kallisto count file contains MAG_ID, length/effective length, counts, and tpm. Initialize empty data structures to receive this data.
dat_kallisto_tpm <- data.frame(bin = character(0))
dat_kallisto_raw <- data.frame(bin = character(0))
dat_kallisto_bin_stats <- data.frame()

# for each individual file, read the contents, extract relevant columns, and use join functions to
# aggregate by column. The MAG lists in each file should be the same, but in case they are not
# we use a fully join.
pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) :eta", total = length(files))
for (i in 1:length(files)) {
  pb$tick()
  SID <- gsub("_kallisto_mapping$", "", gsub("NZ_", "", strsplit(files[i], "/")[[1]][2])) # extract the SID from the filename. The convention may change based on file locations/naming strategies.
  f <- read.table(files[i], sep = "\t", header = TRUE)
  if (i == 1) {
    dat_kallisto_bin_stats <- f %>%
      select(target_id, length, eff_length) %>%
      dplyr::rename(bin = target_id)
  }
  dat_kallisto_tpm <- dat_kallisto_tpm %>%
    full_join(f %>%
                select(target_id, tpm) %>%
                dplyr::rename(bin = target_id, !!SID := tpm), by = "bin")
  dat_kallisto_raw <- dat_kallisto_raw %>%
    full_join(f %>%
                select(target_id, est_counts) %>%
                dplyr::rename(bin = target_id, !!SID := est_counts), by = "bin")
  
}

# some formatting. We're moving away from the "NZ_/NC_" prefix scheme for internal accessions.
dat_kallisto_bin_stats <- dat_kallisto_bin_stats %>%
  mutate(bin = gsub("NZ_", "", bin))

# this is just alignment of SIDs between the MAG abundance datasets and the available metadata
keep <- data.frame(SID = colnames(dat_kallisto_raw)) %>%
  inner_join(map, by = "SID") %>%
  inner_join(map.anthropometry %>% select(-SID), by = c("PID", "study_week"))

# filter the abundance datasets to those samples with matching metadata
dat_kallisto_tpm <- dat_kallisto_tpm %>%
#  rownames_to_column("bin") %>%
  mutate(bin = gsub("^NZ_", "", bin)) %>%
  column_to_rownames("bin") %>%
  select(keep$SID)

dat_kallisto_raw <- dat_kallisto_raw %>%
#  rownames_to_column("bin") %>%
  mutate(bin = gsub("^NZ_", "", bin)) %>%
  column_to_rownames("bin") %>%
  select(keep$SID)

# count data abundance and prevalence filtering ---------------------------

# we often "parameter sweep" to get a feel for the sparsity of our datasets. This code sweeps across reasonable
# values for our presence/absence and prevalence filters and stores the number of MAGs and % total data retained
# it's CRITICAL you consider the samples you intend to analyze (if a subset) when you perform this filtration

# note that here we're filtering based on a TPM-normalized dataset that allows us to more fairly compare MAG abundances
# to each other. Below we'll filter the "raw" (not-normalized) dataset to match the filtered TPM dataset

# set up an empty data frame to receive results
dat_kallisto_tpm_filt.stat <- data.frame(threshold = numeric(0), 
                                         prevalence = numeric(0),
                                         ntaxa = numeric(0),
                                         fraction = numeric(0))

# loop through and filter the abundance dataset
for (threshold in seq(0, 100, 5)) {
  for (prevalence in seq(0, 1, 0.1)) {
    dat_kallisto_tpm.filt <- dat_kallisto_tpm[rowSums(dat_kallisto_tpm > threshold) > (ncol(dat_kallisto_tpm) * prevalence),] %>%
      select(keep$SID)
    dat_kallisto_tpm_filt.stat <- bind_rows(dat_kallisto_tpm_filt.stat,
                                            data.frame(threshold = threshold,
                                                       prevalence = prevalence,
                                                       ntaxa = nrow(dat_kallisto_tpm.filt),
                                                       fraction_zero = sum(dat_kallisto_tpm.filt == 0),
                                                       fraction_data = ifelse(nrow(dat_kallisto_tpm.filt) > 0, sum(dat_kallisto_tpm.filt)/sum(dat_kallisto_tpm), 0)))
  }
}

# plot the grids of the parameter(s) sweep and the results obtained

# fraction data
ggplot(dat_kallisto_tpm_filt.stat, aes(x = prevalence, y = threshold, label = round(fraction_data, 2))) +
  ggtitle("fraction of data remaining") +
  geom_tile(color = "black", fill = NA) +
  geom_text() +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  theme_bw()

# fraction taxa
ggplot(dat_kallisto_tpm_filt.stat, aes(x = prevalence, y = threshold, label = ntaxa)) +
  ggtitle("fraction of taxa remaining") +
  geom_tile(color = "black", fill = NA) +
  geom_text() +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  theme_bw()

# fraction zero
ggplot(dat_kallisto_tpm_filt.stat, aes(x = prevalence, y = threshold, label = round(fraction_zero, 2))) +
  ggtitle("fraction of zeros") +
  geom_tile(color = "black", fill = NA) +
  geom_text() +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  theme_bw()

# select a set of reasonable filtering parameters for your analysis, and apply it to the abundance dataset
# you should have your intended analysis in mind here...

dat_kallisto_tpm.filt <- dat_kallisto_tpm[rowSums(dat_kallisto_tpm > 5) > (ncol(dat_kallisto_tpm) * 0.40),]
sum(dat_kallisto_tpm.filt)/sum(dat_kallisto_tpm) # check that your filter result matches your expectation
nrow(dat_kallisto_tpm.filt) # check that your filter result matches your expectation

# align the "raw" dataset with the filtered dataset
dat_kallisto_raw <- dat_kallisto_raw %>%
  rownames_to_column("MAG") %>%
  filter(MAG %in% rownames(dat_kallisto_tpm.filt)) %>%
  column_to_rownames("MAG")

# gathering some annotation stats about each MAG from prokka output -------

# list the ".txt" summary files from prokka
files <- list.files(path = "annotation_stats", pattern = ".mod2.txt", full.names = TRUE)

# loop through the files, grabbing MAG_ID, coding sequence count, and MAG length
mag_lengths <- data.frame(MAG = character(0),
                          CDS = numeric(0),
                          kb = numeric(0))

for (i in 1:length(files)) {
  mag <- gsub(".txt$", "", gsub("NZ_", "", strsplit(files[i], "/")[[1]][2]))
  mag_prokka_summary <- read.table(files[i], sep = " ", skip = 1, col.names = c("stat", "value"))
  mag_lengths <- mag_lengths %>%
    bind_rows(data.frame(MAG = mag,
                         CDS = mag_prokka_summary[mag_prokka_summary$stat == "CDS:",]$value,
                         kb = mag_prokka_summary[mag_prokka_summary$stat == "bases:",]$value/1000))
}

mag_lengths.anno <- mag_lengths %>%
  left_join(tax %>%
              mutate(MAG = gsub("NZ_", "", MAG)) %>%
              select(MAG, genus, species), by = "MAG") %>%
  arrange(kb)

# write summary to file
write.table(mag_lengths.anno, paste0(datestring, "_MDCF_POC_mag_assembly_lengths.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# abundance data normalization (GETMM) -------------------------------------

# getmm normalization

# calculate RPK
dat_kallisto_rpk <- dat_kallisto_raw %>%
  rownames_to_column("MAG") %>%
  left_join(mag_lengths, by = "MAG") %>%
  column_to_rownames("MAG") %>%
  mutate(across(colnames(dat_kallisto_raw), .fns = function(x) x/kb)) %>%
  select(-kb)

# create an edgeR-friendlyl dataset with a "dummy" grouping factor
rpk.norm <- DGEList(counts = dat_kallisto_rpk,
                              group = c(rep("A",ncol(dat_kallisto_rpk))))

# calculate library size scaling factors via edgeR's TMM approach
rpk.norm <- calcNormFactors(rpk.norm)

# complete the calculation with CPM to achieve GETMM
dat_kallisto_getmm <- cpm(rpk.norm, normalized.lib.sizes = TRUE) %>%
  as.data.frame()

# abundance data normalization (VST) ---------------------------------------

# grab the relevant mapping data
map.sort <- map[match(colnames(dat_kallisto_raw), map$SID),] %>%
  filter(SID %in% colnames(dat_kallisto_raw)) %>% 
  column_to_rownames("SID")
map.sort$study_arm <- factor(map.sort$study_arm, levels = c("X", "Y"))

# prepare a DESeq2-ready data structure
dds <- DESeqDataSetFromMatrix(countData = round(dat_kallisto_raw, 0),
                              colData = map.sort,
                              design = ~ 1)

# estimate library size scaling factors
dds <- estimateSizeFactors(dds, type = "poscounts")

# perform the VST - this is a built-in function in DESeq2
vst <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "local")

# extract the normalized data from the DESeq2 object
dat_kallisto_vst <- assay(vst) %>%
  as.data.frame()

# write normalized datasets to file
write.table(dat_kallisto_raw, paste0(datestring, "_MDCF_POC_MAG_1000_counts_raw_full_set.txt"), quote = FALSE, row.names = TRUE, sep = "\t")
write.table(dat_kallisto_tpm, paste0(datestring, "_MDCF_POC_MAG_1000_counts_tpm_full_set.txt"), quote = FALSE, row.names = TRUE, sep = "\t")
write.table(dat_kallisto_getmm, paste0(datestring, "_MDCF_POC_MAG_1000_counts_getmm_full_set.txt"), quote = FALSE, row.names = TRUE, sep = "\t")
write.table(dat_kallisto_vst, paste0(datestring, "_MDCF_POC_MAG_1000_counts_vst_filt.txt"), quote = FALSE, row.names = TRUE, sep = "\t")

# filter taxonomic data to match filtered datasets, then write reduced set tof ile
tax.filt <- tax %>%
  filter(bin %in% dat_kallisto_bin_stats$bin) %>%
  dplyr::rename(MAG = bin) %>%
  select(MAG, colnames(tax)[!colnames(tax) %in% c("bin", "user_genome")])

write.table(tax.filt, paste0(datestring, "_MDCF_POC_MAG_1000_gtdb_tk_taxonomy.txt"), quote = FALSE, row.names = TRUE, sep = "\t")
