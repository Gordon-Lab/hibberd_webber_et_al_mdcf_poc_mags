library(tidyverse)
library(readxl)
library(ape)
library(cluster)
library(fgsea)
library(edgeR)
library(variancePartition)
library(data.table)
library(naturalsort)
library(scales)

# configure environment
datestring <- format(Sys.time(), "%y%m%d_%H%M")

setwd("~/Documents/Projects/human_studies/MDCF POC/210608_bin_abundance_agg_and_calc_build4/")

# load metadata
map <- read_excel("~/Documents/Projects/human_studies/MDCF POC/200206 16S ASV Analysis/200428_MDCF_POC_fecal_plasma_metadata.xlsx") %>%
  filter(sample_type == "fecal")

map.anthropometry <- read.table("~/Documents/Projects/human_studies/MDCF POC/200423_MDCF_POC_anthropometry.txt", sep = "\t", header = T)

map.quartiles <- read_excel("~/Library/CloudStorage/Box-Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/210810_1323_MDCF_POC_participant_beta_WLZ_quartiles.xlsx", skip = 1) %>%
  select(PID, coef, beta_quartile) %>%
  left_join(map %>% select(PID, study_arm), by = "PID")

map.quartiles %>% 
  group_by(study_arm, beta_quartile) %>%
  summarize(mean(coef))

# load abundance data
dat_kallisto_vst <- read.table("~/Library/CloudStorage/Box-Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/MAGs/210614_0934_MDCF_POC_MAG_1000_counts_vst_filt.txt", sep = "\t", header = TRUE)

# merge metadata and abundance data
dat_kallisto_vst.anthropometry <- dat_kallisto_vst %>%
  rownames_to_column("bin") %>%  
  pivot_longer(-bin, names_to = "SID") %>% 
  pivot_wider(names_from=bin, values_from=value) %>%
  #  rownames_to_column("SID") %>%
  inner_join(map, by = "SID") %>%
  inner_join(map.anthropometry %>% select(PID, study_week, wlz), by = c("PID", "study_week"))

dat_kallisto_vst.anthropometry$study_week <- as.numeric(as.character(dat_kallisto_vst.anthropometry$study_week))

delta_POC <- dat_kallisto_vst.anthropometry %>%
  pivot_longer(c(-PID, -study_week, -SID, -gender, -enrollment_age, -baseline_age, -study_arm, -age_at_sampling, -study_phase_alt, -sample_number, -study_phase, -sample_type, -wlz), names_to = "bin", values_to = "ab") %>%
  pivot_wider(id_cols = c(study_arm, PID, bin), names_from = study_week, values_from = ab)

# fixing variable names
name.fix <- colnames(delta_POC)[!colnames(delta_POC) %in% c("study_arm", "PID", "bin")]
colnames(delta_POC)[!colnames(delta_POC) %in% c("study_arm", "PID", "bin")] <- paste0("wk_", name.fix)

# long form
delta_POC.melt <- delta_POC %>%
  left_join(map.quartiles %>% select(-study_arm), by = "PID") %>%
  pivot_longer(c(-PID, -study_arm, -coef, -beta_quartile, -bin), names_to = "study_week", values_to = "ab")

# loading taxonomy data
tax1 <- read.table("~/Documents/Projects/human_studies/MDCF POC/210115_kallisto_mag_purity/gtdb_checkm_magpurify_clean/gtdbtk.bac120.summary.tsv", sep = "\t", header = TRUE) %>%
  separate(classification, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
  mutate(bin = gsub("^NZ_", "", user_genome))
tax2 <- read.table("~/Documents/Projects/human_studies/MDCF POC/210608_bin_abundance_agg_and_calc_build4/210610_gtdb_tk_33_assignments/gtdbtk.bac120.summary.tsv", sep = "\t", header = TRUE) %>%
  separate(classification, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
  mutate(bin = gsub("^NZ_", "", user_genome))
tax <- bind_rows(tax1, tax2) %>%
  mutate(genus = gsub("g__", "", genus),
         species = gsub("s__", "", species),
         MAG = bin) %>%
  select(-user_genome)

# lists of MAGs
mag_1000 <- read.table("~/Library/CloudStorage/Box-Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/MAGs/210614_0934_MDCF_POC_MAG_1000_counts_tpm_full_set.txt", header = TRUE, sep = "\t") %>%
  rownames_to_column("MAG") %>%
  select(MAG) %>%
  as.vector()

mag_837 <- read.table("~/Library/CloudStorage/Box-Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/MAGs/210614_0934_bin_1000_wlz_res_vst.txt", header = TRUE, sep = "\t") %>%
  select(MAG) %>%
  as.vector()

# results of WLZ-mag relationship testing
res <- read.table("~/Library/CloudStorage/Box-Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/MAGs/210614_0934_bin_1000_wlz_res_vst.txt", header = TRUE, sep = "\t")
sig <- res %>%
  filter(lme_ab_p.adj < 0.05)

# loading mcSEED BPMs
dat <- read_excel("~/Documents/Composition and Presentation/200924 Hibberd Webber et al/versions/230508/from_dmitry/230510_updated_PC_phenotypes.xlsx")

# BPM key
key <- read_excel("../210302_mcSEED_recon_MAGs/mcseed_phenotype_key.xlsx") %>%
  mutate(category = factor(category),
         label = paste(category, name, mcseed_module, sep = "|")) %>%
  arrange(label)

# subset to only those taxa with BPMs
res.anno <- res %>%
  select(MAG, lme_ab_beta, lme_ab_p.adj) %>%
  left_join(tax %>% select(MAG, genus, species), by = "MAG") %>%
  inner_join(dat %>% select(-genus, -species), by = "MAG") %>%
  mutate(sig = ifelse(lme_ab_p.adj <= 0.05, "sig", "ns")) %>%
#  group_by(genus) %>%
#  arrange(sig, .by_group = TRUE)
  arrange(lme_ab_beta)

phenotypes <- colnames(dat)[! colnames(dat) %in% c("species", "genus", "MAG", "Group_Size", "MAG_ID", "beta_WLZ", "q_value")]

# calculations by arm

# loading "raw" unnormalized abundance data
dat_kallisto_raw <- read.table("~/Library/CloudStorage/Box-Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/MAGs/210614_0934_MDCF_POC_MAG_1000_counts_raw_filt.txt", header = TRUE)
metadata <- map %>%
  filter(SID %in% colnames(dat_kallisto_raw),
         study_week < 16) %>%
  inner_join(map.anthropometry %>% select(PID, study_week, wlz, muac, waz, laz), by = c("PID", "study_week")) %>%
  column_to_rownames("SID") %>%
  mutate(study_arm = factor(study_arm, levels = c("Y", "X")))

dat_kallisto_raw <- dat_kallisto_raw %>%
  select(rownames(metadata))

# Standard usage of limma/voom
dge <- DGEList(counts = dat_kallisto_raw,
               samples = metadata)
dge <- calcNormFactors(dge)

# The variable to be tested must be a fixed effect
form <- ~ study_arm + study_week + study_arm:study_week + (1|PID)
#form <- ~ wlz + study_arm + study_week + study_arm:study_week + (1|PID)
#design <- model.matrix(~ study_arm + study_week + study_arm:study_week, data = metadata)

# estimate weights using linear mixed model of dream
vobjDream <- voomWithDreamWeights(counts = dge, 
                                  formula = form, 
                                  data = metadata,
                                  span = 0.8, 
                                  plot = TRUE)

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
fitmm <- dream(exprObj = vobjDream,
               formula = form, 
               data = metadata)

# Examine design matrix
head(fitmm$design, 3)

fitmm.se <- data.frame(fitmm$stdev.unscaled*fitmm$sigma) %>%
  rownames_to_column("MAG")

# Get results of hypothesis test on coefficients of interest
res.study_arm <- topTable(fitmm, coef='study_armX', adjust.method = "BH", n = Inf) %>%
  rownames_to_column("MAG") %>%
  arrange(desc(t)) %>%
  left_join(fitmm.se %>% select(MAG, study_armX) %>% dplyr::rename(se = study_armX), by = "MAG")
res.arm_x_week <- topTable(fitmm, coef='study_armX:study_week', adjust.method = "BH", n = Inf) %>%
  rownames_to_column("MAG") %>%
  arrange(desc(t)) %>%
  left_join(fitmm.se %>% select(MAG, study_armX.study_week) %>% dplyr::rename(se = study_armX.study_week), by = "MAG")
res.study_week <- topTable(fitmm, coef='study_week', adjust.method = "BH", n = Inf) %>%
  rownames_to_column("MAG") %>%
  arrange(desc(t)) %>%
  left_join(fitmm.se %>% select(MAG, study_week) %>% dplyr::rename(se = study_week), by = "MAG")

# enrichment by WLZ association
res.anno.long <- res.anno %>%
  select(-beta_WLZ, -q_value, -sig) %>%
  ungroup() %>%
  select(-lme_ab_beta, -lme_ab_p.adj, -genus, -species) %>%
  pivot_longer(c(-MAG, -MAG_ID), names_to = "mcseed_module", values_to = "phenotype") %>%
  left_join(key, by = "mcseed_module")

# generating the module "sets" for GSEA
mcseed_module_groupings <- list()

for (i in 1:length(phenotypes)) {
  res.anno.long.sub <- res.anno.long %>%
    filter(mcseed_module == phenotypes[i]) %>%
    mutate(phenotype = as.numeric(phenotype))

    mcseed_module_groupings[[paste(phenotypes[i], "1", sep = "_")]] <- res.anno.long.sub %>%
      filter(phenotype == 1) %>%
      pull(MAG)
    mcseed_module_groupings[[paste(phenotypes[i], "0", sep = "_")]] <- res.anno.long.sub %>%
      filter(phenotype == 0) %>%
      pull(MAG)
}

# ranks from WLZ association
ranks <- res.anno$lme_ab_beta
names(ranks) <- res.anno$MAG

# run GSEA
fgsea_res_wlz <- fgsea(mcseed_module_groupings,
                       minSize = 10,
                       maxSize = length(ranks)/2,
                       stats = ranks,
                       eps = 0)

fgsea_res_wlz <- fgsea_res_wlz %>%
  arrange(desc(NES))
fgsea_res_wlz %>%
  filter(padj < 0.05)

# write results
fwrite(fgsea_res_wlz, file = paste0(datestring, "_mcseed_phenotype_by_wlz_fgsea_res.txt"), sep = "\t", sep2 = c("", " ", ""))

# analysis by quartile

# Get results of hypothesis test on coefficients of interest
res.quartile <- read.table("210811_1616_MDCF_POC_MAGs_dream_quartile.txt", sep = "\t", header = TRUE) %>%
  arrange(desc(t))
res.interaction <- read.table("210811_1616_MDCF_POC_MAGs_dream_interaction.txt", sep = "\t", header = TRUE) %>%
  arrange(desc(t))

# checking directions
ranks <- res.quartile$t
names(ranks) <- res.quartile$MAG

fgsea_res_quartile_mcseed <- fgsea(mcseed_module_groupings,
                   minSize = 10,
                   maxSize = length(ranks)/2,
                   stats = ranks,
                   eps = 0)

fgsea_res_quartile_mcseed <- fgsea_res_quartile_mcseed %>%
  arrange(desc(NES))
fgsea_res_quartile_mcseed %>%
  filter(padj < 0.05)

fwrite(fgsea_res_quartile_mcseed, file = paste0(datestring, "_phenotype_by_quartile_MDCF2only_fgsea_res.txt"), sep = "\t", sep2 = c("", " ", ""))

# fgsea of mcseed phenotypes, ranks are interaction t statistic

ranks <- res.interaction$t
names(ranks) <- res.interaction$MAG

# run gsea
fgsea_res_quartile_x_week_mcseed <- fgsea(mcseed_module_groupings,
                   minSize = 10,
                   maxSize = length(ranks)/2,
                   stats = ranks,
                   eps = 0)

fgsea_res_quartile_x_week_mcseed <- fgsea_res_quartile_x_week_mcseed %>%
  arrange(desc(NES))
fgsea_res_quartile_x_week_mcseed %>%
  filter(padj < 0.05)

# write results to file
fwrite(fgsea_res_quartile_x_week_mcseed, file = paste0(datestring, "_phenotype_by_quartile-time_interaction_MDCF2only_fgsea_res.txt"), sep = "\t", sep2 = c("", " ", ""))

# Plotting

sig_phenotypes <- fgsea_res_quartile_mcseed %>%
  select(pathway, padj, NES) %>%
  filter(padj < 0.05) %>%
  mutate(phenotype = gsub("_1|_0", "", pathway)) %>%
  pull(phenotype)
fgsea_plot.quartile <- fgsea_res_quartile_mcseed %>%
  mutate(phenotype = gsub("_1|_0", "", pathway),
         pathway = substr(pathway, nchar(pathway), nchar(pathway))) %>%
  filter(phenotype %in% sig_phenotypes) %>%
  select(pathway, padj, NES, phenotype)

for (i in 1:length(unique(sig_phenotypes))) {
  test2_sub <- fgsea_plot.quartile %>%
    filter(phenotype %in% unique(sig_phenotypes)[i])
  if (nrow(test2_sub) == 2) {
    next()
  } else if ("1" %in% test2_sub$pathway) {
    fgsea_plot.quartile <- fgsea_plot.quartile %>%
      bind_rows(data.frame(pathway = "0",
                           padj = NA,
                           NES = NA,
                           phenotype = unique(sig_phenotypes)[i]))
  } else {
    fgsea_plot.quartile <- fgsea_plot.quartile %>%
      bind_rows(data.frame(pathway = "1",
                           padj = NA,
                           NES = NA,
                           phenotype = unique(sig_phenotypes)[i]))
  }
}

fgsea_plot.quartile <- fgsea_plot.quartile %>%
  left_join(key, by = c("phenotype" = "mcseed_module")) %>%
  group_by(phenotype) %>%
  mutate(order = sum(abs(NES), na.rm = TRUE)) %>%
  mutate(sig = ifelse(padj < 0.05 & !is.na(padj), "sig", "ns"),
         label = paste(category, name, phenotype, sep = "|"))
fgsea_plot.quartile$label <- factor(fgsea_plot.quartile$label, levels = rev(key$label[key$label %in% fgsea_plot.quartile$label]))

# by interaction

sig_phenotypes <- fgsea_res_quartile_x_week_mcseed %>%
  select(pathway, padj, NES) %>%
  filter(padj < 0.05) %>%
  mutate(phenotype = gsub("_1|_0", "", pathway)) %>%
  pull(phenotype)
fgsea_plot.quartile_x_week <- fgsea_res_quartile_x_week_mcseed %>%
  mutate(phenotype = gsub("_1|_0", "", pathway),
         pathway = substr(pathway, nchar(pathway), nchar(pathway))) %>%
  filter(phenotype %in% sig_phenotypes) %>%
  select(pathway, padj, NES, phenotype)

for (i in 1:length(unique(sig_phenotypes))) {
  test2_sub <- fgsea_plot.quartile_x_week %>%
    filter(phenotype %in% unique(sig_phenotypes)[i])
  if (nrow(test2_sub) == 2) {
    next()
  } else if ("1" %in% test2_sub$pathway) {
    fgsea_plot.quartile_x_week <- fgsea_plot.quartile_x_week %>%
      bind_rows(data.frame(pathway = "0",
                           padj = NA,
                           NES = NA,
                           phenotype = unique(sig_phenotypes)[i]))
  } else {
    fgsea_plot.quartile_x_week <- fgsea_plot.quartile_x_week %>%
      bind_rows(data.frame(pathway = "1",
                           padj = NA,
                           NES = NA,
                           phenotype = unique(sig_phenotypes)[i]))
  }
}

fgsea_plot.quartile_x_week <- fgsea_plot.quartile_x_week %>%
  left_join(key, by = c("phenotype" = "mcseed_module")) %>%
  group_by(phenotype) %>%
  mutate(order = sum(abs(NES), na.rm = TRUE)) %>%
  mutate(sig = ifelse(padj < 0.05 & !is.na(padj), "sig", "ns"),
         label = paste(category, name, phenotype, sep = "|"))
fgsea_plot.quartile_x_week$label <- factor(fgsea_plot.quartile_x_week$label, levels = rev(key$label[key$label %in% fgsea_plot.quartile_x_week$label]))

# by wlz_association

sig_phenotypes <- fgsea_res_wlz %>%
  select(pathway, padj, NES) %>%
  filter(padj < 0.05) %>%
  mutate(phenotype = gsub("_1|_0", "", pathway)) %>%
  pull(phenotype)
fgsea_plot.wlz <- fgsea_res_wlz %>%
  mutate(phenotype = gsub("_1|_0", "", pathway),
         pathway = substr(pathway, nchar(pathway), nchar(pathway))) %>%
  filter(phenotype %in% sig_phenotypes) %>%
  select(pathway, padj, NES, phenotype)

for (i in 1:length(unique(sig_phenotypes))) {
  test2_sub <- fgsea_plot.wlz %>%
    filter(phenotype %in% unique(sig_phenotypes)[i])
  if (nrow(test2_sub) == 2) {
    next()
  } else if ("1" %in% test2_sub$pathway) {
    fgsea_plot.wlz <- fgsea_plot.wlz %>%
      bind_rows(data.frame(pathway = "0",
                           padj = NA,
                           NES = NA,
                           phenotype = unique(sig_phenotypes)[i]))
  } else {
    fgsea_plot.wlz <- fgsea_plot.wlz %>%
      bind_rows(data.frame(pathway = "1",
                           padj = NA,
                           NES = NA,
                           phenotype = unique(sig_phenotypes)[i]))
  }
}

fgsea_plot.wlz <- fgsea_plot.wlz %>%
  left_join(key, by = c("phenotype" = "mcseed_module")) %>%
  group_by(phenotype) %>%
  mutate(order = sum(abs(NES), na.rm = TRUE)) %>%
  mutate(sig = ifelse(padj < 0.05 & !is.na(padj), "sig", "ns"),
         label = paste(category, name, phenotype, sep = "|"))
fgsea_plot.wlz$label <- factor(fgsea_plot.wlz$label, levels = rev(key$label[key$label %in% fgsea_plot.wlz$label]))

# by arm

ranks <- res.study_arm$t
names(ranks) <- res.study_arm$MAG

fgsea_res_arm_mcseed <- fgsea(mcseed_module_groupings,
                              minSize = 10,
                              maxSize = length(ranks)/2,
                              stats = ranks,
                              eps = 0)

fgsea_res_arm_mcseed <- fgsea_res_arm_mcseed %>%
  arrange(desc(NES))
fgsea_res_arm_mcseed %>%
  filter(padj < 0.05)

fwrite(fgsea_res_arm_mcseed, file = paste0(datestring, "_phenotype_by_arm_fgsea_res.txt"), sep = "\t", sep2 = c("", " ", ""))

sig_phenotypes <- fgsea_res_arm_mcseed %>%
  select(pathway, padj, NES) %>%
  filter(padj < 0.05) %>%
  mutate(phenotype = gsub("_1|_0", "", pathway)) %>%
  pull(phenotype)
fgsea_plot.arm <- fgsea_res_arm_mcseed %>%
  mutate(phenotype = gsub("_1|_0", "", pathway),
         pathway = substr(pathway, nchar(pathway), nchar(pathway))) %>%
  filter(phenotype %in% sig_phenotypes) %>%
  select(pathway, padj, NES, phenotype)

for (i in 1:length(unique(sig_phenotypes))) {
  test2_sub <- fgsea_plot.arm %>%
    filter(phenotype %in% unique(sig_phenotypes)[i])
  if (nrow(test2_sub) == 2) {
    next()
  } else if ("1" %in% test2_sub$pathway) {
    fgsea_plot.arm <- fgsea_plot.arm %>%
      bind_rows(data.frame(pathway = "0",
                           padj = NA,
                           NES = NA,
                           phenotype = unique(sig_phenotypes)[i]))
  } else {
    fgsea_plot.arm <- fgsea_plot.arm %>%
      bind_rows(data.frame(pathway = "1",
                           padj = NA,
                           NES = NA,
                           phenotype = unique(sig_phenotypes)[i]))
  }
}

fgsea_plot.arm <- fgsea_plot.arm %>%
  left_join(key, by = c("phenotype" = "mcseed_module")) %>%
  group_by(phenotype) %>%
  mutate(order = sum(abs(NES), na.rm = TRUE)) %>%
  mutate(sig = ifelse(padj < 0.05 & !is.na(padj), "sig", "ns"),
         label = paste(category, name, phenotype, sep = "|"))
fgsea_plot.arm$label <- factor(fgsea_plot.arm$label, levels = rev(key$label[key$label %in% fgsea_plot.arm$label]))

# by interaction

ranks <- res.interaction$t
names(ranks) <- res.interaction$MAG

fgsea_res_arm_x_week_mcseed <- fgsea(mcseed_module_groupings,
                                      minSize = 10,
                                      maxSize = length(ranks)/2,
                                      stats = ranks,
                                      eps = 0)

fgsea_res_arm_x_week_mcseed <- fgsea_res_arm_x_week_mcseed %>%
  arrange(desc(NES))
fgsea_res_arm_x_week_mcseed %>%
  filter(padj < 0.05)

fgsea_res_arm_x_week_mcseed %>%
  filter(size < 20)

fwrite(fgsea_res_arm_x_week_mcseed, file = paste0(datestring, "_phenotype_by_interaction_fgsea_res.txt"), sep = "\t", sep2 = c("", " ", ""))

sig_phenotypes <- fgsea_res_arm_x_week_mcseed %>%
  select(pathway, padj, NES) %>%
  filter(padj < 0.05) %>%
  mutate(phenotype = gsub("_1|_0", "", pathway)) %>%
  pull(phenotype)
fgsea_plot.arm_x_week <- fgsea_res_arm_x_week_mcseed %>%
  mutate(phenotype = gsub("_1|_0", "", pathway),
         pathway = substr(pathway, nchar(pathway), nchar(pathway))) %>%
  filter(phenotype %in% sig_phenotypes) %>%
  select(pathway, padj, NES, phenotype)

for (i in 1:length(unique(sig_phenotypes))) {
  test2_sub <- fgsea_plot.arm_x_week %>%
    filter(phenotype %in% unique(sig_phenotypes)[i])
  if (nrow(test2_sub) == 2) {
    next()
  } else if ("1" %in% test2_sub$pathway) {
    fgsea_plot.arm_x_week <- fgsea_plot.arm_x_week %>%
      bind_rows(data.frame(pathway = "0",
                           padj = NA,
                           NES = NA,
                           phenotype = unique(sig_phenotypes)[i]))
  } else {
    fgsea_plot.arm_x_week <- fgsea_plot.arm_x_week %>%
      bind_rows(data.frame(pathway = "1",
                           padj = NA,
                           NES = NA,
                           phenotype = unique(sig_phenotypes)[i]))
  }
}

fgsea_plot.arm_x_week <- fgsea_plot.arm_x_week %>%
  left_join(key, by = c("phenotype" = "mcseed_module")) %>%
  group_by(phenotype) %>%
  mutate(order = sum(abs(NES), na.rm = TRUE)) %>%
  mutate(sig = ifelse(padj < 0.05 & !is.na(padj), "sig", "ns"),
         label = paste(category, name, phenotype, sep = "|"))
fgsea_plot.arm_x_week$label <- factor(fgsea_plot.arm_x_week$label, levels = rev(key$label[key$label %in% fgsea_plot.arm_x_week$label]))

fgsea_plot.agg <- fgsea_plot.wlz %>%
  ungroup() %>%
  select(pathway, label, NES, sig) %>%
  mutate(stat_group = "wlz") %>%
  bind_rows(fgsea_plot.arm %>%
              ungroup() %>%
              select(pathway, label, NES, sig) %>%
              mutate(stat_group = "arm")) %>%
  bind_rows(fgsea_plot.arm_x_week %>%
              ungroup() %>%
              select(pathway, label, NES, sig) %>%
              mutate(stat_group = "arm_x_week")) %>%
  bind_rows(fgsea_plot.quartile %>%
              ungroup() %>%
              select(pathway, label, NES, sig) %>%
              mutate(stat_group = "quartile")) %>%
  bind_rows(fgsea_plot.quartile_x_week %>%
              ungroup() %>%
              select(pathway, label, NES, sig) %>%
              mutate(stat_group = "quartile_x_week")) %>%
  mutate(stat_group = factor(stat_group, levels = c("wlz", "arm", "arm_x_week", "quartile", "quartile_x_week")))

#  filter(pathway != "0")

fgsea_plot.agg$label <- factor(fgsea_plot.agg$label, levels = rev(key$label[key$label %in% fgsea_plot.agg$label]))

for (i in 1:length(unique(fgsea_plot.agg$label))) {
  test2_sub <- fgsea_plot.agg %>%
    filter(label %in% unique(fgsea_plot.agg$label)[i])
  if (nrow(test2_sub) == 10) {
    next()
  } else {
    for (g in c("wlz", "arm", "arm_x_week", "quartile", "quartile_x_week")) {
      test3_sub <- test2_sub %>%
        filter(stat_group == g)
      if (nrow(test3_sub) == 2) {
        next()
      } else {
        if (! "1" %in% test3_sub$label) {
          fgsea_plot.agg <- fgsea_plot.agg %>%
            bind_rows(data.frame(stat_group = g,
                                 pathway = "1",
                                 #                                 padj = NA,
                                 NES = NA,
                                 label = unique(test2_sub$label)))
        }
        if (! "0" %in% test3_sub$label) {
          fgsea_plot.agg <- fgsea_plot.agg %>%
            bind_rows(data.frame(stat_group = g,
                                 pathway = "0",
                                 #                               padj = NA,
                                 NES = NA,
                                 label = unique(test2_sub$label)))
        }
      }
    }
  }
}

# aggregate plots across multiple analyses

pdf(paste0(datestring, "_phenotype_agg_grey_bright.pdf"))
ggplot(fgsea_plot.agg, aes(y = label, x = pathway, fill = NES)) +
  geom_tile(color = "black") +
  geom_tile(aes(height = 1, width = 1), data = fgsea_plot.agg %>% filter(sig == "sig"), color = "black") +
  coord_equal() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "grey95") +
  theme_classic() +
  facet_grid(. ~ stat_group)
dev.off()

summary <- fgsea_plot.agg %>%
  mutate(sign = ifelse(is.na(NES), "ns", ifelse(NES > 0, "pos", "neg")),
         sign_sig = paste0(sign, "_", sig)) %>%
  pivot_wider(id_cols = c(label, pathway), names_from = stat_group, values_from = sign_sig) %>%
  mutate(wlz_arm_con = ifelse(arm == wlz & arm != "ns_NA", TRUE, FALSE),
         wlz_arm_x_week_con = ifelse(arm_x_week == wlz & arm_x_week != "ns_NA", TRUE, FALSE),
         wlz_quartile_con = ifelse(quartile == wlz & quartile != "ns_NA", TRUE, FALSE),
         wlz_quartile_x_week_con = ifelse(quartile_x_week == wlz & quartile_x_week != "ns_NA", TRUE, FALSE)) %>%
  filter(wlz_arm_con == TRUE | wlz_arm_x_week_con == TRUE | wlz_quartile_con == TRUE | wlz_quartile_x_week_con == TRUE)

fgsea_plot.agg.filt <- fgsea_plot.agg %>%
  filter(label %in% summary$label) %>%
  mutate(stat_group = factor(stat_group, levels = c("wlz", "arm", "arm_x_week", "quartile", "quartile_x_week")),
         pathway = factor(as.character(pathway), levels = c("1", "0")))

pdf(paste0(datestring, "_phenotype_agg_grey_bright.pdf"))
ggplot(fgsea_plot.agg.filt, aes(y = label, x = pathway, fill = NES)) +
  geom_tile(color = "black") +
  geom_tile(aes(height = 1, width = 1), data = fgsea_plot.agg.filt %>% filter(sig == "sig"), color = "black") +
  coord_equal() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "grey95") +
  theme_classic() +
  facet_grid(. ~ stat_group)
dev.off()

pdf(paste0(datestring, "_phenotype_agg_grey_dark.pdf"))
ggplot(fgsea_plot.agg.filt, aes(y = label, x = pathway, fill = NES)) +
  geom_tile(color = "black") +
  geom_tile(aes(height = 1, width = 1), data = fgsea_plot.agg.filt %>% filter(sig == "sig"), color = "black") +
  coord_equal() +
  scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), na.value = "grey95") +
  theme_classic() +
  facet_grid(. ~ stat_group)
dev.off()
