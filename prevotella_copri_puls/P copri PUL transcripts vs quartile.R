# configure environment

library(tidyverse)
library(DESeq2)
library(readxl)
library(naturalsort)
library(edgeR)
library(variancePartition)
library(naturalsort)

datestring <- format(Sys.time(), "%y%m%d_%H%M")
study <- "MDCF_POC_PM"

setwd("~/Documents/Projects/human_studies/MDCF POC/220502 PULs vs MDCF-2 WLZ MAGs/")

# metadata
map <- read_excel("~/Documents/Projects/human_studies/MDCF POC/200206 16S ASV Analysis/200428_MDCF_POC_fecal_plasma_metadata.xlsx") %>%
  filter(sample_type == "fecal")
map.anthropometry <- read.table("~/Documents/Projects/human_studies/MDCF POC/200423_MDCF_POC_anthropometry.txt", sep = "\t", header = T) %>%
  mutate(PID = factor(PID))
map.quartiles <- read_excel("~/Library/CloudStorage/Box-Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/210810_1323_MDCF_POC_participant_beta_WLZ_quartiles.xlsx", skip = 1) %>%
  select(PID, coef, beta_quartile) %>%
  left_join(map %>% select(PID, study_arm), by = "PID") %>%
  distinct()

# check quartile directions
map.quartiles %>% 
  group_by(study_arm, beta_quartile) %>%
  summarize(mean(coef))

load("~/Documents/Projects/human_studies/MDCF POC/220502 PULs vs MDCF-2 WLZ MAGs/210621_MDCF_POC_transcript_counts_raw_filtered_wTranscriptZeroing.Rda")

dat_rna <- counts_ind_flt %>%
  rownames_to_column("locus_tag")
rm(counts_ind_flt)

dat_rna <- dat_rna %>%
  mutate(locus_base = gsub("_.*", "", locus_tag))

# lookup for locus tag bases for each MAG
bin_tag_key <- read.table("~/Documents/Projects/human_studies/MDCF POC/210608_bin_abundance_agg_and_calc_build4/210609_bin_tag_key.txt", col.names = c("MAG", "locus_base")) %>%
  mutate(MAG = gsub("^NZ_", "", MAG),
         locus_base = gsub("META", "MBIN", locus_base))
mag_info <- read_excel("~/Library/CloudStorage/Box-Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/code for deposition/hibberd_webber_et_al_mdcf_poc_mags/cc_svd/CodeFigure4AB/210924_table_s1_mag_info.xlsx") %>%
  left_join(bin_tag_key, by = c("mag_name" = "MAG")) %>%
  mutate(mag_id = gsub("_", "", mag_id))

target_mags <- c("MAGBg0018", "MAGBg0019")
pc_mags <- c("MAGBg0018", "MAGBg0019", "MAGBg0008", "MAGBg0003", "MAGBg0004", "MAGBg0002", "MAGBg0048", "MAGBg0035", "MAGBg0015", "MAGBg0101", "MAGBg0138")
target_puls <- read_excel("220512_SM_pc_puls.xlsx")

pc_tags <- mag_info %>%
#  filter(mag_id %in% pc_mags) %>%
  filter(mag_id %in% "MAGBg0018") %>%
  pull(locus_base)

# subset to P. copri MAGs
dat_rna.sub <- dat_rna %>%
  filter(locus_base %in% pc_tags)

# subset to P. copri of interest and their component PUL transcripts
dat_rna.sub.puls <- dat_rna.sub %>%
  inner_join(target_puls, by = c("locus_tag")) %>%
  pivot_longer(cols = c(-locus_tag, -locus_base, -pul, -MAG_ID, -consensus_pul, -conservation, -CAZyme, -gene, -protein_product), names_to = "SID", values_to = "count") %>%
#  mutate(SID = gsub("_[P|Q]$", "", SID)) %>%
  select(-locus_base)

all_participants <- unique(map$PID)

# merge metadata and PUL transcript data
dat_rna.sub.puls.anno <- dat_rna.sub.puls %>%
  inner_join(map, by = "SID") %>%
  inner_join(map.anthropometry %>% select(PID, study_week, wlz), by = c("PID", "study_week"))

# summarize
dat_rna.sub.puls.anno.summary <- dat_rna.sub.puls.anno %>%
  mutate(phase_week = paste(study_phase, study_week, sep = "_")) %>%
  group_by(MAG_ID, locus_tag, pul, consensus_pul, conservation, study_arm, phase_week) %>%
  summarize(count_mean = mean(count),
            count_sd = sd(count))

# getting a feel for directions
dat_rna.sub.puls.anno.summary <- dat_rna.sub.puls.anno.summary %>%
  mutate(count_mean = ifelse(count_mean == 0, NA, count_mean),
         phase_week = factor(phase_week, levels = c("baseline_0", "treatment_4", "treatment_12")))

# some plotting for PULs of interest
ggplot(dat_rna.sub.puls.anno.summary %>% filter(MAG_ID %in% target_mags & consensus_pul %in% c("PUL_7", "PUL_8")), aes(x = phase_week, y = locus_tag, fill = count_mean)) +
  geom_tile() +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred") +
  facet_grid(consensus_pul ~ study_arm, scales = "free_y")

ggplot(dat_rna.sub.puls.anno.summary, aes(x = phase_week, y = locus_tag, fill = count_mean)) +
  geom_tile() +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred") +
  facet_grid(. ~ study_arm)

# data normalization

#'full' vst
dat_raw <- dat_rna %>%
  column_to_rownames("locus_tag") %>%
  select(-locus_base)

#vst on subset MAGs assigned to P. copri
dat_raw <- dat_rna.sub %>%
  column_to_rownames("locus_tag") %>%
  select(-locus_base)
colnames(dat_raw) <- gsub("_[P|Q]$", "", colnames(dat_raw))

# check for zero columns
zero_cols <- colnames(dat_raw)[colSums(dat_raw) == 0]
dat_raw <- dat_raw %>%
  select(-all_of(zero_cols))

# align metadata
map.sort <- map[match(colnames(dat_raw), map$SID),] %>%
  filter(SID %in% colnames(dat_raw)) %>% 
  column_to_rownames("SID")
map.sort$study_arm <- factor(map.sort$study_arm, levels = c("X", "Y"))

dds <- DESeqDataSetFromMatrix(countData = round(dat_raw, 0),
                              colData = map.sort,
                              design = ~ 1)

dds <- estimateSizeFactors(dds, type = "poscounts")
vst <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "local")
dat_vst <- assay(vst) %>%
  as.data.frame()

# VST-normalized data for PUL transcripts
dat_vst.puls <- dat_vst %>%
  rownames_to_column("locus_tag") %>%
  mutate(locus_base = gsub("_.*", "", locus_tag)) %>%
  inner_join(target_puls, by = c("locus_tag"))

# add metadata
dat_vst.puls.anno <- dat_vst.puls %>%
  pivot_longer(cols = c(-locus_tag, -locus_base, -pul, -MAG_ID, -consensus_pul, -conservation, -CAZyme, -gene, -protein_product), names_to = "SID", values_to = "count") %>%
  inner_join(map, by = "SID") %>%
  inner_join(map.anthropometry %>% select(PID, study_week, wlz), by = c("PID", "study_week")) %>%
  left_join(map.quartiles %>% select(PID, beta_quartile), by = "PID")
  
# summarize by arm
dat_vst.puls.anno.summary.arm <- dat_vst.puls.anno %>%
  mutate(phase_week = paste(study_phase, study_week, sep = "_")) %>%
  group_by(MAG_ID, locus_tag, pul, consensus_pul, conservation, study_arm, phase_week) %>%
  summarize(count_mean = mean(count),
            count_sd = sd(count)) %>%
  mutate(phase_week = factor(phase_week, levels = c("baseline_0", "treatment_4", "treatment_12")),
         locus_tag = factor(locus_tag, levels = naturalsort(unique(locus_tag))))

# plotting
ggplot(dat_vst.puls.anno.summary %>% filter(MAG_ID %in% target_mags & consensus_pul %in% c("PUL_7", "PUL_8")), aes(x = phase_week, y = locus_tag, fill = count_mean)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred") +
  facet_grid(consensus_pul ~ study_arm, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5))

# summarize by quartile
dat_vst.puls.anno.summary.quartile <- dat_vst.puls.anno %>%
  mutate(phase_week = paste(study_phase, study_week, sep = "_")) %>%
  group_by(MAG_ID, locus_tag, pul, consensus_pul, conservation, study_arm, beta_quartile, phase_week) %>%
  summarize(count_mean = mean(count),
            count_sd = sd(count)) %>%
  mutate(phase_week = factor(phase_week, levels = c("baseline_0", "treatment_4", "treatment_12")),
         locus_tag = factor(locus_tag, levels = naturalsort(unique(locus_tag))))

# plot by quartile, for puls of interest, upper vs lower quartile
ggplot(dat_vst.puls.anno.summary.quartile %>% filter(MAG_ID %in% target_mags & consensus_pul %in% c("PUL_7", "PUL_8") & beta_quartile %in% c("Q1", "Q4") & study_arm == "X"), aes(x = phase_week, y = locus_tag, fill = count_mean)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred") +
  facet_grid(consensus_pul ~ beta_quartile, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5))

dat_vst.puls.anno.summary.quartile <- dat_vst.puls.anno %>%
  mutate(phase_week = paste(study_phase, study_week, sep = "_")) %>%
  select(MAG_ID, PID, locus_tag, pul, consensus_pul, conservation, study_arm, beta_quartile, phase_week, count) %>%
  distinct() %>%
  pivot_wider(id_cols = c(PID, MAG_ID, locus_tag, pul, consensus_pul, conservation, study_arm, beta_quartile), names_from = phase_week, values_from = count) %>%
  mutate(delta = treatment_12 - baseline_0) %>%
#  summarize(mean_delta = mean(delta),
#            sd_delta = sd(delta)) %>%
  select(-baseline_0, -treatment_4, -treatment_12) %>%
  pivot_longer(delta, values_to = "delta") %>%
  mutate(locus_tag = factor(locus_tag, levels = naturalsort(unique(locus_tag))))

ggplot(dat_vst.puls.anno.summary.quartile %>% filter(MAG_ID %in% target_mags & consensus_pul %in% c("PUL_7", "PUL_8", "PUL_4", "PUL_17b") & beta_quartile %in% c("Q1", "Q4") & study_arm == "X"), aes(x = delta, y = locus_tag, fill = MAG_ID)) +
#  geom_point() +
  geom_vline(xintercept = 0) +
  geom_boxplot() +
  facet_grid(consensus_pul ~ beta_quartile) +
  theme_classic()

pdf(paste(datestring, study, "vst_norm_pul_transcription_all_participants.pdf", sep = "_"), height = 8.5, width = 11)
ggplot(dat_vst.puls.anno.summary, aes(x = phase_week, y = locus_tag, fill = count_mean)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred") +
  facet_grid(consensus_pul ~ study_arm, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5))
dev.off()

# dream - for extended data plotting

metadata.treatment <- map %>%
  left_join(map.quartiles %>% select(PID, beta_quartile), by = "PID") %>%
  filter(SID %in% colnames(dat_raw) & study_phase %in% c("baseline", "treatment")) %>% 
#  filter(beta_quartile %in% c("Q1", "Q4") & study_arm == "X") %>%
  mutate(study_arm = factor(study_arm, levels = c("Y", "X"))) %>%
  column_to_rownames("SID")

dat_raw.treatment <- dat_raw %>%
  select(rownames(metadata.treatment))

# Standard usage of limma/voom
dge.treatment <- DGEList(counts = dat_raw.treatment,
                         samples = metadata.treatment)
dge.treatment <- calcNormFactors(dge.treatment)

# The variable to be tested must be a fixed effect
form.treatment <- ~  beta_quartile*study_week + (1|PID)
form.treatment <- ~  study_arm*study_week + (1|PID)
#design <- model.matrix(~ study_arm + study_week + study_arm:study_week, data = metadata)

# estimate weights using linear mixed model of dream
dream_obj.treatment <- voomWithDreamWeights(counts = dge.treatment, 
                                            formula = form.treatment, 
                                            data = metadata.treatment,
                                            span = 0.8, 
                                            plot = TRUE,
                                            save.plot = FALSE)

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
dream_fit.treatment <- dream(exprObj = dream_obj.treatment,
                             formula = form.treatment, 
                             data = metadata.treatment)

#dream_fit.treatment <- eBayes(dream_fit.treatment)

# differences by interaction
res.int <- topTable(dream_fit.treatment, coef='study_armX:study_week', adjust.method = "BH", n = Inf) %>%
  rownames_to_column("locus_tag") %>%
  arrange(z.std)

res.int <- topTable(dream_fit.treatment, coef='beta_quartileQ4:study_week', adjust.method = "BH", n = Inf) %>%
  rownames_to_column("locus_tag") %>%
  arrange(z.std)

# for pul transcripts of interest
res.int.puls <- res.int %>%
  filter(locus_tag %in% target_puls$locus_tag)

# add metadata
res.int.puls.anno <- res.int.puls %>%
  left_join(target_puls, by = "locus_tag") %>%
  filter(! is.na(consensus_pul)) %>%
  distinct() %>%
  arrange(pul, locus_tag) %>%
  mutate(locus_tag = factor(locus_tag, levels = locus_tag))

# write results to file
write.table(res.int.puls.anno, "220818_quartile_x_week_mag18_puls.txt", sep = "\t", row.names = FALSE)
write.table(res.int.puls.anno, "220818_diet_x_week_mag18_puls.txt", sep = "\t", row.names = FALSE)

# plotting
ggplot(res.int.puls.anno, aes(x = z.std, y = locus_tag)) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(x = 0, xend = z.std, y = locus_tag, yend = locus_tag)) +
  geom_point() +
  facet_grid(consensus_pul ~ ., scales = "free_y") +
  theme_classic()
  
# normalizing for mag abundance

dat_dna.raw <- read.table("~/Library/CloudStorage/Box-Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/MAGs/210614_0934_MDCF_POC_MAG_1000_counts_raw_filt.txt", header = TRUE)

map.sort.dna <- map %>%
  filter(SID %in% colnames(dat_dna.raw),
         study_week < 16) %>%
  left_join(map.quartiles %>% select(PID, beta_quartile), by = "PID") %>%
  column_to_rownames("SID") %>%
  mutate(study_arm = factor(study_arm, levels = c("Y", "X")))

dat_dna.raw <- dat_dna.raw %>%
  select(rownames(map.sort.dna))

dds.dna <- DESeqDataSetFromMatrix(countData = round(dat_dna.raw, 0),
                              colData = map.sort.dna,
                              design = ~ 1)

dds.dna <- estimateSizeFactors(dds.dna, type = "poscounts")
vst.dna <- varianceStabilizingTransformation(dds.dna, blind = TRUE, fitType = "local")
dat_dna.vst <- assay(vst.dna) %>%
  as.data.frame()

save.image("220513_1210_pcopri_focus_pul_transcription.RData")

# plotting

dat_dream.arm <- read.table("220818_diet_x_week_mag19_puls.txt", sep = "\t", header = TRUE)
dat_dream.quartile <- read.table("220818_quartile_x_week_mag19_puls.txt", sep = "\t", header = TRUE)

dat_dream.plot <- dat_dream.arm %>%
  mutate(analysis = "arm_x_time") %>%
  bind_rows(dat_dream.quartile %>%
              mutate(analysis = "quartile_x_time")) %>%
  mutate(label = paste(locus_tag, CAZyme, protein_product, sep = "|")) %>%
  arrange(locus_tag) %>%
  mutate(label = factor(label, levels = unique(label)))

ggplot(dat_dream.plot, aes(x = z.std, y = locus_tag)) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(x = 0, xend = z.std, y = locus_tag, yend = locus_tag)) +
  geom_point() +
  facet_grid(consensus_pul ~ analysis, scales = "free_y") +
  theme_classic()

ggplot(dat_dream.plot %>% filter(consensus_pul %in% c("PUL_7", "PUL_8", "PUL_4", "PUL_17b", "PUL_17a")), aes(x = logFC, y = locus_tag)) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(x = 0, xend = logFC, y = locus_tag, yend = locus_tag)) +
  geom_point() +
  facet_grid(consensus_pul ~ analysis, scales = "free_y") +
  theme_classic()

ggplot(dat_dream.plot %>% filter(consensus_pul %in% c("PUL_7", "PUL_8", "PUL_4", "PUL_17b", "PUL_17a")), aes(x = logFC, y = label)) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(x = 0, xend = logFC, y = label, yend = label)) +
  geom_point() +
  facet_grid(consensus_pul ~ analysis, scales = "free") +
  theme_classic()


