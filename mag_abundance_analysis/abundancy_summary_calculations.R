# configure environment
library(tidyverse)
library(readxl)

setwd("/Users/hibberdm/Documents/Projects/human_studies/MDCF POC/210608_bin_abundance_agg_and_calc_build4")

datestring <- format(Sys.time(), "%y%m%d_%H%M")

# load VST and TPM-normalized abundance data
dat_kallisto_vst <- read.table("~/Box/MDCF_POC_Microbiome_Paper/MAGs/210614_0934_MDCF_POC_MAG_1000_counts_vst_filt.txt", header = TRUE, sep = "\t")
dat_kallisto_tpm <- read.table("~/Box/MDCF_POC_Microbiome_Paper/MAGs/210614_0934_MDCF_POC_MAG_1000_counts_tpm_full_set.txt", header = TRUE, sep = "\t") %>%
  rownames_to_column("MAG") %>%
  filter(MAG %in% rownames(dat_kallisto_vst)) %>%
  select(c(MAG, colnames(dat_kallisto_vst))) %>%
  column_to_rownames("MAG")

# load metadata
map <- read_excel("~/Documents/Projects/human_studies/MDCF POC/200206 16S ASV Analysis/200428_MDCF_POC_fecal_plasma_metadata.xlsx") %>%
  filter(sample_type == "fecal")

map.anthropometry <- read.table("~/Documents/Projects/human_studies/MDCF POC/200423_MDCF_POC_anthropometry.txt", sep = "\t", header = T) %>%
  mutate(PID = factor(PID))

# load taxonomy
tax1 <- read.table("~/Documents/Projects/human_studies/MDCF POC/210115_kallisto_mag_purity/gtdb_checkm_magpurify_clean/gtdbtk.bac120.summary.tsv", sep = "\t", header = TRUE) %>%
  separate(classification, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
  mutate(MAG = gsub("^NZ_", "", user_genome))
tax2 <- read.table("~/Documents/Projects/human_studies/MDCF POC/210608_bin_abundance_agg_and_calc_build4/210610_gtdb_tk_33_assignments/gtdbtk.bac120.summary.tsv", sep = "\t", header = TRUE) %>%
  separate(classification, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
  mutate(MAG = gsub("^NZ_", "", user_genome))
tax <- bind_rows(tax1, tax2) %>%
  mutate(genus = gsub("g__", "", genus),
         species = gsub("s__", "", species))

# merge metadata with abundance data
dat_kallisto_vst.anthropometry <- dat_kallisto_vst %>%
  rownames_to_column("bin") %>%  
  pivot_longer(-bin, names_to = "SID") %>% 
  pivot_wider(names_from=bin, values_from=value) %>%
  #  rownames_to_column("SID") %>%
  merge(map, by = "SID") %>%
  merge(map.anthropometry %>% select(PID, study_week, wlz), by = c("PID", "study_week"))

# df formatting
delta_POC <- dat_kallisto_vst.anthropometry %>%
  pivot_longer(c(-PID, -study_week, -SID, -gender, -enrollment_age, -baseline_age, -study_arm, -age_at_sampling, -study_phase_alt, -sample_number, -study_phase, -sample_type, -wlz), names_to = "bin", values_to = "ab") %>%
  pivot_wider(id_cols = c(study_arm, PID, bin), names_from = study_week, values_from = ab)

# fix variable names
name.fix <- colnames(delta_POC)[!colnames(delta_POC) %in% c("study_arm", "PID", "bin")]
colnames(delta_POC)[!colnames(delta_POC) %in% c("study_arm", "PID", "bin")] <- paste0("wk_", name.fix)

# long form
delta_POC.melt <- delta_POC %>%
  left_join(beta_table_v2 %>% select(-study_arm), by = "PID") %>%
  pivot_longer(c(-PID, -study_arm, -coef, -delta_wlz, -beta_quartile, -bin), names_to = "study_week", values_to = "ab")

# calculate summary statistics
delta_POC.melt.stat <- delta_POC.melt %>%
  group_by(bin, study_arm, study_week) %>%
  mutate(study_week = as.numeric(gsub("wk_", "", study_week))) %>%
  summarize(mean = mean(ab, na.rm = TRUE), 
            sd = sd(ab, na.rm = TRUE), 
            se = sd(ab, na.rm = TRUE) / sqrt(length(ab)),
            prev_0 = sum(ab > 0)) %>%
  mutate(mean = round(mean, 1),
         sd = round(sd, 1),
         se = round(se, 1),
         msd = paste(mean, "±", sd, sep = " "),
         mse = paste(mean, "±", se, sep = " "))

# delta_POC.melt.stat %>%
#   pivot_wider(id_cols = c(bin, study_arm), names_from = study_week, values_from = msd) %>%
#   arrange(study_arm, bin) %>%
#   left_join(tax %>% select(bin, genus, species), by = "bin") %>%
#   write.table(paste0(datestring, "_MDCF_POC_MAG_abundance_summary_msd.txt"), row.names = FALSE, sep = "\t")
# 
# delta_POC.melt.stat %>%
#   pivot_wider(id_cols = c(bin, study_arm), names_from = study_week, values_from = mse) %>%
#   arrange(study_arm, bin) %>%
#   left_join(tax %>% select(bin, genus, species), by = "bin") %>%
#   write.table(paste0(datestring, "_MDCF_POC_MAG_abundance_summary_mse.txt"), row.names = FALSE, sep = "\t")
# 
# delta_POC.melt.stat %>%
#   pivot_wider(id_cols = c(bin, study_arm), names_from = study_week, values_from = prev_0) %>%
#   arrange(study_arm, bin) %>%
#   left_join(tax %>% select(bin, genus, species), by = "bin") %>%
#   write.table(paste0(datestring, "_MDCF_POC_MAG_abundance_summary_prev_0.txt"), row.names = FALSE, sep = "\t")

#write.table(wk_0_data, paste0(datestring, "_MDCF_POC_MAG_abundance_summary_week_0.txt"), row.names = FALSE, sep = "\t")  

# TPM data for prevalence reporting

# merge abundance data and metadata
dat_kallisto_tpm.anthropometry <- dat_kallisto_tpm %>%
  rownames_to_column("bin") %>%  
  pivot_longer(-bin, names_to = "SID") %>% 
  pivot_wider(names_from=bin, values_from=value) %>%
  #  rownames_to_column("SID") %>%
  merge(map %>% select(SID, PID, study_week, study_arm), by = "SID")
#  merge(map.anthropometry %>% select(PID, study_week, wlz), by = c("PID", "study_week"))

# df formatting
delta_POC_tpm <- dat_kallisto_tpm.anthropometry %>%
  select(-SID) %>%
  pivot_longer(c(-PID, -study_week, -study_arm), names_to = "bin", values_to = "ab") %>%
  pivot_wider(id_cols = c(study_arm, PID, bin), names_from = study_week, values_from = ab)

# fix variable names
name.fix <- colnames(delta_POC_tpm)[!colnames(delta_POC_tpm) %in% c("study_arm", "PID", "bin")]
colnames(delta_POC_tpm)[!colnames(delta_POC_tpm) %in% c("study_arm", "PID", "bin")] <- paste0("wk_", name.fix)

# long form
delta_POC_tpm.melt <- delta_POC_tpm %>%
#  left_join(beta_table_v2 %>% select(-study_arm), by = "PID") %>%
  pivot_longer(c(-PID, -study_arm, -bin), names_to = "study_week", values_to = "ab")

# summary stats
delta_POC_tpm.melt.stat <- delta_POC_tpm.melt %>%
  group_by(bin, study_arm, study_week) %>%
  mutate(study_week = as.numeric(gsub("wk_", "", study_week))) %>%
  summarize(mean = mean(ab, na.rm = TRUE), 
            sd = sd(ab, na.rm = TRUE), 
            se = sd(ab, na.rm = TRUE) / sqrt(length(ab)),
            prev_5 = sum(ab[!is.na(ab)] > 5)) %>%
  mutate(mean = round(mean, 1),
         sd = round(sd, 1),
         se = round(se, 1),
         msd = paste(mean, "±", sd, sep = " "),
         mse = paste(mean, "±", se, sep = " "))

# prevalence table by MAG
prev_table <- delta_POC_tpm.melt.stat %>%
  pivot_wider(id_cols = c(bin, study_arm), names_from = study_week, values_from = prev_5) %>%
  arrange(study_arm, bin) %>%
  left_join(tax %>% select(MAG, genus, species), by = c("bin" = "MAG"))

#%>%
#  write.table(paste0(datestring, "_MDCF_POC_MAG_abundance_summary_tpm_prev5.txt"), row.names = FALSE, sep = "\t")

# working on prevalence statistics
delta_POC_tpm.melt.stat %>%
  pivot_wider(id_cols = c(bin, study_arm), names_from = study_week, values_from = msd) %>%
  arrange(study_arm, bin) %>%
  left_join(tax %>% select(MAG, genus, species), by = c("bin" = "MAG")) %>%
  write.table(paste0(datestring, "_MDCF_POC_MAG_abundance_summary_tpm_msd.txt"), row.names = FALSE, sep = "\t")

delta_POC.melt.stat <- delta_POC_tpm.melt %>%
  group_by(bin, study_week) %>%
  mutate(study_week = as.numeric(gsub("wk_", "", study_week))) %>%
  summarize(mean = mean(ab, na.rm = TRUE), 
            sd = sd(ab, na.rm = TRUE), 
            se = sd(ab, na.rm = TRUE) / sqrt(length(ab)),
            prev_0 = sum(ab > 0)) %>%
  mutate(mean = round(mean, 1),
         sd = round(sd, 1),
         se = round(se, 1),
         msd = paste(mean, "±", sd, sep = " "),
         mse = paste(mean, "±", se, sep = " "))

wk_0_data <- delta_POC.melt.stat %>%
  filter(study_week == 0) %>%
  arrange(desc(mean))
wk_0_data$MAG_ID = paste0("MAG_", formatC(seq(1, 837, 1), width = 4, format = "d", flag = "0"))

## calculating average MAGs count per individual per week per arm

# prevalence table
prev_table <- delta_POC_tpm.melt %>%
  group_by(study_arm, PID, study_week) %>%
  summarize(mag_gt_5 = sum(ab > 100, na.rm = TRUE)) %>%
  pivot_wider(id_cols = c(study_arm, PID), names_from = study_week, values_from = mag_gt_5) %>%
  group_by(study_arm) %>%
  summarize(across(starts_with("wk_"), .fns = list(mean, sd)))

# plotting
ggplot(prev_table, aes(x = study_week, y = mag_gt_5)) +
  geom_boxplot() +
  facet_grid(. ~ study_arm)

# long form
delta_POC_tpm.melt %>%
  group_by(study_arm, PID, study_week) %>%
  summarize(mag_gt_5 = sum(ab > 5, na.rm = TRUE)) %>%
#  pivot_wider(id_cols = c(study_arm, PID), names_from = study_week, values_from = mag_gt_5) %>%
  group_by(study_arm) %>%
  summarize(across(mag_gt_5, .fns = list(mean, sd)))

# summarizing MAG assembly stats

# load quast results
mag_assembly_stats <- read.table("quast_results/latest/transposed_report.tsv", skip = 1, sep = "\t")
colnames(mag_assembly_stats) <- c("MAG", "contigs_gt0", "contigs_gt1000", "contigs_gt5000", "contigs_gt10000", "contigs_gt25000", "contigs_gt50000",
                         "length_gt0", "length_gt1000", "length_gt5000", "length_gt10000", "length_gt25000", "length_gt50000",
                         "contigs", "longest_contig", "length", "gc", "N50", "N75", "L50", "L75", "N_per_100kbp")
mag_assembly_stats <- mag_assembly_stats %>%
  mutate(MAG = gsub("NZ_", "", MAG))

# load checkm results
checkm_res <- read.table("210614_MDCF_POC_MAG_1000_checkm_res.txt", header = FALSE, sep = "\t", skip = 1)
colnames(checkm_res) <- c("bin",
                          "marker_lineage",
                          "n_genomes",
                          "n_markers",
                          "n_marker_sets",
                          "n_0",
                          "n_1",
                          "n_2",
                          "n_3",
                          "n_4",
                          "n_5plus",
                          "completeness",
                          "contamination",
                          "strain_heterogeneity")

checkm_res <- checkm_res %>%
  mutate(MAG = gsub("NZ_", "", bin)) %>%
  select(-bin)

# load prokka summary stats
files <- list.files(path = "annotation_stats", pattern = ".mod2.txt", full.names = TRUE)
mag_annotation_stats <- data.frame()

for (i in 1:length(files)) {
  mag <- gsub(".txt$", "", gsub("NZ_", "", strsplit(files[i], "/")[[1]][2]))
  mag_prokka_summary <- read.table(files[i], sep = " ", skip = 1, col.names = c("stat", "value")) %>%
    mutate(stat = gsub(":$", "", stat)) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    mutate(length = bases/1000,
           MAG = mag)

  mag_annotation_stats <- mag_annotation_stats %>%
    bind_rows(mag_prokka_summary)
    
}

# merge all summary tables
table_S2 <- wk_0_data %>%
  dplyr::rename(MAG = bin) %>%
  left_join(mag_annotation_stats, by = "MAG") %>%
  left_join(mag_assembly_stats %>% select(MAG, contigs, length, gc, N50), by = "MAG") %>%
  left_join(checkm_res %>% select(MAG, completeness, contamination, strain_heterogeneity), by = "MAG") %>%
  left_join(tax, by = "MAG") %>%
  left_join(wlz_asv_evidence, by = "MAG") %>%
  mutate(kingdom = gsub("d__", "", kingdom),
         phylum = gsub("p__", "", phylum),
         class = gsub("c__", "", class),
         order = gsub("o__", "", order),
         family = gsub("f__", "", family))

# write result to file
write.table(table_S2, paste0(datestring, "_MDCF_POC_microbiome_table_S2.txt"), row.names = FALSE, sep = "\t")  
