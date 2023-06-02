# configure environment
library(fgsea)
library(tidyverse)
library(readxl)
library(progress)
library(ggrepel)
library(progress)
library(DESeq2)
library(edgeR)
library(variancePartition)
library(gridExtra)

setwd("~/Documents/Projects/human_studies/MDCF POC/211014 Prevotella focus based on PULs/")

datestring <- format(Sys.time(), "%y%m%d_%H%M")
study <- "MDCF_POC"
experiment <- "wlz-MAGs"

cbbPalette <- c("#999999", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "darkred", "darkblue", "darkgrey", "white")

# stable functions
source("functions.R")

# load metadata
map <- read_excel("~/Documents/Projects/human_studies/MDCF POC/200206 16S ASV Analysis/200428_MDCF_POC_fecal_plasma_metadata.xlsx") %>%
  filter(sample_type == "fecal")

map.anthropometry <- read.table("~/Documents/Projects/human_studies/MDCF POC/200423_MDCF_POC_anthropometry.txt", sep = "\t", header = T)

map.quartiles <- read_excel("~/Box/MDCF_POC_Microbiome_Paper/210810_1323_MDCF_POC_participant_beta_WLZ_quartiles.xlsx", skip = 1) %>%
  select(PID, coef, beta_quartile) %>%
  left_join(map %>% select(PID, study_arm), by = "PID")

# confirm which quartile is "upper"
map.quartiles %>% 
  group_by(study_arm, beta_quartile) %>%
  summarize(mean(coef))

# load taxonomy data
tax <- read.table("../210608_bin_abundance_agg_and_calc_build4/210614_1509_MDCF_POC_MAG_1000_gtdb_tk_taxonomy.txt", sep = "\t", header = TRUE)

# load raw and VST-normalized abundance data
dat_kallisto_raw <- read.table("~/Box/MDCF_POC_Microbiome_Paper/MAGs/210614_0934_MDCF_POC_MAG_1000_counts_raw_filt.txt", sep = "\t", header = TRUE)
dat_kallisto_vst <- read.table("~/Box/MDCF_POC_Microbiome_Paper/MAGs/210614_0934_MDCF_POC_MAG_1000_counts_vst_filt.txt", sep = "\t", header = TRUE)

# merging metadata and abundance data
dat_kallisto_vst.anthropometry <- dat_kallisto_vst %>%
  rownames_to_column("bin") %>%  
  pivot_longer(-bin, names_to = "SID") %>% 
  pivot_wider(names_from=bin, values_from=value) %>%
  #  rownames_to_column("SID") %>%
  inner_join(map, by = "SID") %>%
  inner_join(map.anthropometry %>% select(PID, study_week, wlz), by = c("PID", "study_week"))

dat_kallisto_vst.anthropometry$study_week <- as.numeric(as.character(dat_kallisto_vst.anthropometry$study_week))

# Checking direction

target_mags <- read.table("pco_targets.txt", sep = "\t", col.names = c("short_name", "mag_name"))
#%>%
#  filter(short_name %in% c("MAG_Bg_0018", "MAG_Bg_0019", "MAG_Bg_0048"))

# orienting dataset
delta_POC <- dat_kallisto_vst.anthropometry %>%
  pivot_longer(c(-PID, -study_week, -SID, -gender, -enrollment_age, -baseline_age, -study_arm, -age_at_sampling, -study_phase_alt, -sample_number, -study_phase, -sample_type, -wlz), names_to = "bin", values_to = "ab") %>%
  pivot_wider(id_cols = c(study_arm, PID, bin), names_from = study_week, values_from = ab)

# fixing names
name.fix <- colnames(delta_POC)[!colnames(delta_POC) %in% c("study_arm", "PID", "bin")]
colnames(delta_POC)[!colnames(delta_POC) %in% c("study_arm", "PID", "bin")] <- paste0("wk_", name.fix)

# long form
delta_POC.melt <- delta_POC %>%
  left_join(map.quartiles %>% select(-study_arm), by = "PID") %>%
  pivot_longer(c(-PID, -study_arm, -coef, -beta_quartile, -bin), names_to = "study_week", values_to = "ab")

# positive by quartile
dat.plot_pc <- delta_POC.melt %>%
  filter(bin %in% target_mags$mag_name) %>%
#         study_arm == "X",
#         beta_quartile %in% c("Q1", "Q4")) %>%
  left_join(target_mags, by = c("bin" = "mag_name")) %>%
  mutate(study_week = factor(study_week, levels = c("wk_0", "wk_2", "wk_4", "wk_8", "wk_12", "wk_16")))

ggplot(dat.plot_pc, aes(x = study_week, y = ab, color = study_arm)) +
#  geom_point() +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~ short_name)

ggplot(dat.plot_pc, aes(x = study_week, y = ab, color = study_arm, group = study_arm)) +
  geom_smooth(method = "lm") +
  theme_bw() +
  facet_wrap(~ short_name)

ggplot(dat.plot_pc, aes(x = study_week, y = ab, color = study_arm, group = study_arm)) +
  geom_smooth(method = "loess") +
  theme_bw() +
  facet_wrap(~ short_name)

ggplot(dat.plot_pc %>% filter(study_week != "wk_16"), aes(x = study_week, y = ab, color = study_arm, group = study_arm)) +
  geom_smooth(method = "lm") +
  theme_bw() +
  facet_wrap(~ short_name)

##

dat.plot_pc2 <- delta_POC.melt %>%
  filter(bin %in% target_mags$mag_name,
#         study_arm == "Y",
         beta_quartile %in% c("Q1")) %>%
  left_join(target_mags, by = c("bin" = "mag_name")) %>%
  mutate(study_week = factor(study_week, levels = c("wk_0", "wk_2", "wk_4", "wk_8", "wk_12", "wk_16")))

ggplot(dat.plot_pc2, aes(x = study_week, y = ab, color = study_arm, group = study_arm)) +
  geom_smooth(method = "loess") +
  theme_classic() +
  facet_wrap(~ short_name, scales = "free_y")

# 
# load wlz mag info from statistical analysis of WLZ~MAG abundance associations
res <- read.table("~/Box/MDCF_POC_Microbiome_Paper/MAGs/210614_0934_bin_1000_wlz_res_vst.txt", header = TRUE, sep = "\t")
sig <- res %>%
  filter(lme_ab_p.adj < 0.05)
res.target <- res %>% 
  inner_join(target_mags, by = c("MAG" = "mag_name")) %>%
  arrange(desc(lme_ab_beta)) %>%
  mutate(short_name = gsub("_", "", short_name),
         short_name = factor(short_name, levels = short_name),
         sig = ifelse(lme_ab_p.adj < 0.05, TRUE, FALSE))

# pul conservation analysis

# load pul conservation data
pul_conservation <- read_excel("~/Documents/Projects/human_studies/MDCF POC/220106 Prevotella PUL CAZyme update/220106_pul_conservation_for_clustering.xlsx", sheet = "coded")

# remove non-conserved PULs
pul_conservation <- pul_conservation %>%
  filter(!consensus_pul %in% c("3", "16", "17a", "17b"))

# calculate a euclidean distance between MAGs/isolates, based on PUL conservation
pul_dist <- pul_conservation %>% 
  select(-substrate) %>% 
  pivot_longer(-consensus_pul, names_to = "mag") %>% 
  pivot_wider(names_from=consensus_pul, values_from=value) %>%
  column_to_rownames("mag") %>%
  dist(method = "euclidean")

# hierarchical clustering for organization
pul_clust <- hclust(pul_dist, method = "complete")

# convert to data frame
pul_dist.df <- as.data.frame(as.matrix(pul_dist)) %>%
  rownames_to_column("mag_1") %>%
  pivot_longer(-mag_1, names_to = "mag_2") %>%
  filter(mag_1 == "MAGBg0019" & mag_2 != "PS131S11") %>%
  left_join(res.target %>%
              select(lme_ab_beta, short_name, lme_ab_p.adj),
            by = c("mag_2" = "short_name"))
  
# correlating PUL conservation with WLZ-association
cor.test(pul_dist.df$value, pul_dist.df$lme_ab_beta)

# plot result
ggplot(pul_dist.df, aes(x = value, y = lme_ab_beta, label = mag_2)) +
  geom_point() +
#  geom_line() +
#  geom_label_repel(nudge_x = 0.5, nudge_y = 0.001) +
  theme_classic()

# plot with labels
ggplot(pul_dist.df, aes(x = value, y = lme_ab_beta, label = mag_2)) +
  geom_smooth(method = "lm") +
  geom_point() +
  geom_label_repel(data = pul_dist.df %>% filter(mag_2 %in% c("MAGBg0018", "MAGBg0019")), nudge_x = 0.5, nudge_y = 0.001) +
  theme_classic()
