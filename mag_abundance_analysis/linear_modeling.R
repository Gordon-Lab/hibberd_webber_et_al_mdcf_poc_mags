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
library(naturalsort)
library(corrplot)

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

# read in taxonomic classifications for filtered MAGs from GTDB-tk.
tax <- read.table("210614_1509_MDCF_POC_MAG_1000_gtdb_tk_taxonomy.txt", sep = "\t", header = TRUE)

# read in quartile assignments
map.quartiles <- read_excel("~/Library/CloudStorage/Box-Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/210810_1323_MDCF_POC_participant_beta_WLZ_quartiles.xlsx", skip = 1) %>%
  select(PID, coef, beta_quartile) %>%
  left_join(map %>% select(PID, study_arm), by = "PID")

# check/orient to the designations for "upper", "lower", etc.
map.quartiles %>% 
  group_by(study_arm, beta_quartile) %>%
  summarize(mean(coef))

# raw and vst-normalized abundance data
dat_kallisto_raw <- read.table("~/Library/CloudStorage/Box-Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/MAGs/210614_0934_MDCF_POC_MAG_1000_counts_raw_filt.txt", sep = "\t", header = TRUE)
dat_kallisto_vst <- read.table("~/Library/CloudStorage/Box-Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/MAGs/210614_0934_MDCF_POC_MAG_1000_counts_vst_filt.txt", sep = "\t", header = TRUE)

# modeling MAG abundances versus anthropometry ----------------------------

# adding metadata/anthropometry data to the abundance dataset dataframe
# now is the time to select the subset of samples for your intended analysis if you haven't already done so at the filtering step above.
dat_kallisto_vst.anthropometry <- dat_kallisto_vst %>%
  rownames_to_column("bin") %>%  
  pivot_longer(-bin, names_to = "SID") %>% 
  pivot_wider(names_from=bin, values_from=value) %>%
  #  rownames_to_column("SID") %>%
  merge(map, by = "SID") %>%
  merge(map.anthropometry %>% select(PID, study_week, wlz), by = c("PID", "study_week"))
#%>%
#  filter(study_week < 16) 

dat_kallisto_vst.anthropometry$study_week <- as.numeric(as.character(dat_kallisto_vst.anthropometry$study_week))

# optional - check the individual MAG abundance data for normality, zero fraction, etc.

# zinf <- data.frame(zero_count = rowSums(dat_kallisto_vst.filt == 0)) %>%
#   arrange(zero_count) %>%
#   rownames_to_column("bin")
# 
# for (i in seq(nrow(zinf), nrow(zinf)-10, -1)) {
#   p <- ggplot(dat_kallisto_vst.anthropometry, aes(x = .data[[zinf$bin[i]]], y = wlz)) +
#     geom_point()
#   print(p)
# }
# 
# zinf$sw_p.val <- NA
# 
# for (i in 1:nrow(zinf)) {
#   zinf$sw_p.val[i] <- shapiro.test(dat_kallisto_vst.anthropometry[[zinf$bin[i]]])$p.value
# }

# define functions that define the intended model. The approach makes it easier to apply these functions and store results
model.f <- function(x) { 
  mag_abundance <- as.numeric(x)
  #  lmerTest::lmer(wlz ~ mag_abundance + study_week + mag_abundance:study_week + (1|PID), data = dat_kallisto_vst.anthropometry)
  lmerTest::lmer(wlz ~ mag_abundance + study_week + (1|PID), data = dat_kallisto_vst.anthropometry) 
}

# we use anova to determine the significance of terms.
anova.f <- function(x) { 
  stats::anova(x, ddf="Kenward-Roger") 
}

# list the MAGs for analysis
mags <- rownames(dat_kallisto_vst)

# run the models and anova analyses. store the results as objects in two lists.
models.f <- list()
anovas.f <- list()
pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) :eta", total = length(mags))
for (i in 1:length(mags)) {
  pb$tick()
  models.f[[mags[i]]] <- model.f(dat_kallisto_vst.anthropometry[[mags[i]]])
  anovas.f[[mags[i]]] <- anova.f(models.f[[mags[i]]])
}

# define terms of interest
ab_term <- "mag_abundance"
rate_term <- "mag_abundance:study_week"

# loop through the models/anovas lists and extract the metrics related to those terms of interest.

lme_res_kallisto.full <- data.frame()

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) :eta", total = length(models.f))
for (i in 1:length(mags)) {
  pb$tick()
  ab_portion <- data.frame(MAG = mags[i],
                           lme_ab_beta = get_coef(models.f[[mags[i]]], ab_term),
                           lme_ab_beta_std = get_std_coef(models.f[[mags[i]]], ab_term, dat_kallisto_vst.anthropometry[[mags[i]]], dat_kallisto_vst.anthropometry$wlz),
                           lme_ab_p.val = get_anova_p(anovas.f[[mags[i]]], ab_term),
                           lme_ab_ci = get_confint(models.f[[mags[i]]], ab_term),
                           lme_ab_se = get_model_se(models.f[[mags[i]]], ab_term),
                           lme_rate_beta = get_coef(models.f[[mags[i]]], rate_term),
                           lme_rate_beta_std = get_std_coef(models.f[[mags[i]]], rate_term, dat_kallisto_vst.anthropometry[[mags[i]]], dat_kallisto_vst.anthropometry$wlz),
                           lme_rate_p.val = get_anova_p(anovas.f[[mags[i]]], rate_term),
                           lme_rate_ci = get_confint(models.f[[mags[i]]], rate_term),
                           lme_rate_se = get_model_se(models.f[[mags[i]]], rate_term))
  
  lme_res_kallisto.full <- lme_res_kallisto.full %>%
    bind_rows(ab_portion)
  
}

# adjust p values across all MAGs
lme_res_kallisto.full <- lme_res_kallisto.full %>%
  mutate(lme_ab_p.adj = p.adjust(lme_ab_p.val, method = "fdr"),
         lme_rate_p.adj = p.adjust(lme_rate_p.val, method = "fdr"))

# how many significant interactions were observed for terms of interest?
sum(lme_res_kallisto.full$lme_ab_p.adj < 0.05)
sum(lme_res_kallisto.full$lme_rate_p.adj < 0.05)

# filter the results to only those deemed 'significant' (this may take different forms)
sig.kallisto_vst <- lme_res_kallisto.full %>%
  #  filter(lme_ab_p.adj <= 0.05 | lme_pa_p.adj <= 0.05) %>%
  filter(lme_ab_p.adj <= 0.05) %>%
  #  filter(lme_ab_beta > 0.1 | lme_ab_beta < -0.1) %>%
  #  slice_max(order_by = abs(lme_ab_beta), prop = 0.25) %>%
  #  mutate(MAG = paste0("NZ_", MAG)) %>%
  left_join(tax %>% select(bin, genus, species), by = c("MAG" = "bin")) %>%
  #  distinct() %>%
  mutate(sign = ifelse(lme_ab_beta > 0, "pos", "neg")) %>%
  arrange(lme_ab_beta)

pos_wlz_mags <- sig.kallisto_vst %>%
  filter(sign == "pos") %>%
  pull(MAG)

neg_wlz_mags <- sig.kallisto_vst %>%
  filter(sign == "neg") %>%
  pull(MAG)


# write results to file
write.table(sig.kallisto_vst, file = paste0(datestring, "_bin_1000_wlz_res_sig_vst_treatment_only.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(lme_res_kallisto.full, file = paste0(datestring, "_bin_1000_wlz_res_vst_treatment_only.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

# checking model fits

pdf("residuals.pdf")
for (i in 1:length(mags)) {
  print(i)
  resid <- data.frame(fitted = fitted(models.f[[mags[i]]]),
                      residual = residuals(models.f[[mags[i]]], type = "pearson", scaled = TRUE))
  p1 <- ggplot(resid, aes(x = fitted, y = residual)) +
    geom_hline(yintercept = 0) +
    geom_point(color = "blue") +
    ggtitle(mags[i]) +
    theme_bw()
  p2 <- ggplot(resid, aes(sample = residual)) +
    #  geom_point(color = "blue") +
    ggtitle(mags[i]) +
    xlab("Theoretical Quantiles") +
    ylab("Sample Quantiles") +
    theme_bw() +
    stat_qq() +
    stat_qq_line()
  grid.arrange(p1, p2, nrow = 2)
}
dev.off()

# save checkpoint for MAG vs anthropometry analysis -----------------------

save.image(paste0(datestring, "_wlz_MAGs_vst_1000_post.RData"))

# running PCA on the MAG abundance profiles -------------------------------

# PCA via svd (missing values are an issue here)
pca <- prcomp(dat_kallisto_vst, center = TRUE, scale = TRUE)

# Other options for how to run PCA
#stats::princomp #spectral decomposition
#stats::prcomp #svd
#FactoMineR::PCA #svd

# via eigendecomposition of the covariance matrix (more robust to missing values)
# pca <- eigen(cov(dat_kallisto_vst))
# pca <- eigen(cov(dat_kallisto_vst %>%
#                    rownames_to_column("MAG") %>%
#                    filter(MAG %in% c(pos_wlz_mags, neg_wlz_mags)) %>%
#                    column_to_rownames("MAG")))
# pca <- eigen(cov(dat_kallisto_vst %>%
#                    rownames_to_column("MAG") %>%
#                    filter(MAG %in% pos_wlz_mags) %>%
#                    column_to_rownames("MAG")))
# pca <- eigen(cov(dat_kallisto_vst %>%
#                    rownames_to_column("MAG") %>%
#                    filter(MAG %in% neg_wlz_mags) %>%
#                    column_to_rownames("MAG")))

# extract the variance explained by each PC
var_expl <- data.frame(label = paste("PC", seq(1, length(pca$values)), sep = ""), var_expl = pca$values/sum(pca$values))
var_expl$label <- naturalfactor(var_expl$label)

# focus on the first three - this is somewhat arbitrary but usually useful
pc1_var <- var_expl$var_expl[1]
pc2_var <- var_expl$var_expl[2]
pc3_var <- var_expl$var_expl[3]
variance_summary <- paste("PC1: ", pc1_var, "PC2: ", pc2_var, "PC3: ", pc3_var, sep = "\t")

# extract the coordinates of each sample along each PC
pca_df <- as.data.frame(pca$vectors)
colnames(pca_df) <- paste("PC", seq(1, ncol(pca_df)), sep = "")

# annotate the coordinates with relevant metadata
pca_df$SID <- colnames(dat_kallisto_vst)
pca_df.anno <- pca_df %>%
  inner_join(map, by = "SID") %>%
  inner_join(map.anthropometry %>% select(PID, study_week, wlz), by = c("PID", "study_week")) %>%
  left_join(map.quartiles %>% distinct() %>% select(PID, beta_quartile), by = "PID")

# scree plot - how many axes capture meaningful variance  
p <- ggplot(var_expl[1:10,], aes(x = label, y = var_expl)) +
  geom_bar(stat = "identity") +
  theme_bw()
plot(p)

# export data
write.table(pca_df, paste0(datestring, "_", study, "_", experiment, "_PCA_res.txt"), row.names = TRUE, col.names = TRUE, sep = "\t", quote = F)
write.table(variance_summary, paste0(datestring, "_", study, "_", experiment, "_PCA_var.txt"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = F)

# prepare for plotting. calculate centroid locations for groups of samples along each PC
pca_df.anno.cent <- pca_df.anno %>%
#  filter(study_arm != "RUSF") %>%
  group_by(study_arm, study_phase_alt) %>%
  summarize_at(vars(starts_with("PC")), mean)

# plot samples
ggplot(pca_df.anno, aes(x = PC1, y = PC2, color = study_phase_alt, shape = study_arm)) +
  geom_point(size = 2) +
  scale_color_manual(values = cbbPalette) +
  theme_bw() +
  xlab(paste0("PC1 (", round(pc1_var*100, digits = 2), "% variance)")) +
  ylab(paste0("PC2 (", round(pc2_var*100, digits = 2), "% variance)"))

# plot centroids
ggplot(pca_df.anno.cent, aes(x = PC1, y = PC2, color = study_phase_alt, shape = study_arm)) +
  geom_point(size = 4) +
  geom_path(aes(group = study_arm)) +
  scale_color_manual(values = cbbPalette) +
  theme_bw() +
  xlab(paste0("PC1 (", round(pc1_var*100, digits = 2), "% variance)")) +
  ylab(paste0("PC2 (", round(pc2_var*100, digits = 2), "% variance)"))
#  guides(color = FALSE)

# group and aggregate by quartile
pca_df.anno.cent2 <- pca_df.anno %>%
  #  filter(study_arm != "RUSF") %>%
  group_by(study_arm, beta_quartile, study_phase_alt) %>%
  summarize_at(vars(starts_with("PC")), mean)

# plot centroids
ggplot(pca_df.anno.cent2, aes(x = PC1, y = PC2, color = study_phase_alt, shape = study_arm)) +
  geom_point(size = 4) +
  geom_path(aes(group = study_arm)) +
  scale_color_manual(values = cbbPalette) +
  theme_bw() +
  xlab(paste0("PC1 (", round(pc1_var*100, digits = 2), "% variance)")) +
  ylab(paste0("PC2 (", round(pc2_var*100, digits = 2), "% variance)")) +
  facet_grid(study_arm ~ beta_quartile)

# save checkpoint for PCA -------------------------------------------------

save.image(paste0(datestring, "_PCA.RData"))

# MAG abundance comparisons at baseline (by group) -----------------------------------

# using dream on "raw" abundance data here

# extract the relevant subset of metadata
metadata.baseline <- map %>%
  filter(SID %in% colnames(dat_kallisto_raw),
         study_week == 0) %>%
  mutate(study_arm = factor(study_arm, levels = c("Y", "X"))) %>%
  column_to_rownames("SID")

# align the raw dataset to the selected metadata
dat_kallisto_raw.baseline <- dat_kallisto_raw %>%
  select(rownames(metadata.baseline))

# prepare a limma/voom/edgeR/dream-compatible data structure
# edgeR doesn't handle mixed effects - dream does but the overall approach is different
dge.baseline <- DGEList(dat_kallisto_raw.baseline)

# calculate library size normalization factors
dge.baseline <- calcNormFactors(dge.baseline)

# estimate precision weights (mean-variance relationships) for each MAG. This approach builds on the 
# 'traditional' limma/voom approach to allow mixed effects. 
vobjDream.baseline <- voomWithDreamWeights(dge.baseline, 
                                           ~ study_arm, 
                                           metadata.baseline)

# Fit the models for each MAG, then perform hypothesis tests as defined in the contrast matrix (default or custom)
# degrees of freedom are estimated using the Satterthwaite approximation by default (Kenward-Roger is more accurate, but slower - recommended for < 10 samples)
fitmm.baseline <- dream(vobjDream.baseline,
                        ~ study_arm, 
                        metadata.baseline)

head(fitmm.baseline$design, 3)

# calculate empirial Bayes moderation of standard error -based statistics and log-odds of DE
fitmm.baseline <- eBayes(fitmm.baseline)

# extract results by study arm
res.baseline <- topTable(fitmm.baseline, coef='study_armX', adjust.method = "BH", n = Inf) %>%
  arrange(desc(t)) %>%
  rownames_to_column("MAG") %>%
  left_join(tax %>% select(MAG, genus, species), by = "MAG")

res.baseline %>%
  filter(adj.P.Val < 0.05)

# MAG abundance comparisons at baseline (by quartile) -----------------------------------

# relevant metadata
metadata.baseline.quartile <- map %>%
  filter(SID %in% colnames(dat_kallisto_raw),
         study_week == 0,
         study_arm == "X") %>%
  left_join(map.quartiles %>% distinct() %>% select(PID, beta_quartile), by = "PID") %>%
  mutate(beta_quartile = factor(beta_quartile, levels = c("Q1", "Q2", "Q3", "Q4"))) %>%
  column_to_rownames("SID")

# align abundance data
dat_kallisto_raw.baseline.quartile <- dat_kallisto_raw %>%
  select(rownames(metadata.baseline.quartile))

# prepare a limma/voom/edgeR/dream-compatible data structure
# edgeR doesn't handle mixed effects - dream does but the overall approach is different
dge.baseline.quartile <- DGEList(dat_kallisto_raw.baseline.quartile)

# calculate library size normalization factors
dge.baseline.quartile <- calcNormFactors(dge.baseline.quartile)

# estimate precision weights (mean-variance relationships) for each MAG. This approach builds on the 
# 'traditional' limma/voom approach to allow mixed effects. 
vobjDream.baseline.quartile <- voomWithDreamWeights(dge.baseline.quartile, 
                                                    ~ beta_quartile, 
                                                    metadata.baseline.quartile)

# Fit the models for each MAG, then perform hypothesis tests as defined in the contrast matrix (default or custom)
# degrees of freedom are estimated using the Satterthwaite approximation by default (Kenward-Roger is more accurate, but slower - recommended for < 10 samples)
fitmm.baseline.quartile <- dream(vobjDream.baseline.quartile,
                                 ~ beta_quartile, 
                                 metadata.baseline.quartile)

# calculate empirial Bayes moderation of standard error -based statistics and log-odds of DE
fitmm.baseline.quartile <- eBayes(fitmm.baseline.quartile)

head(fitmm.baseline.quartile$design, 3)

# extract results by quartiles
res.baseline.Q4Q1 <- topTable(fitmm.baseline.quartile, coef='beta_quartileQ4', adjust.method = "BH", n = Inf) %>%
  arrange(desc(t)) %>%
  rownames_to_column("MAG") %>%
  left_join(tax %>% select(MAG, genus, species), by = "MAG")

res.baseline.Q4Q1 %>%
  filter(adj.P.Val < 0.05)

# modeling MAG abundances as a function of MDCF-2 versus RUSF trea --------

# optional: if starting from this point in the code, load the data
#dat_kallisto_raw <- read.table("210614_0934_MDCF_POC_MAG_1000_counts_raw_filt.txt", header = TRUE)

# subset the metadata to the intended analysis
metadata <- map %>%
  filter(SID %in% colnames(dat_kallisto_raw),
         study_week < 16) %>%
  column_to_rownames("SID") %>%
  mutate(study_arm = factor(study_arm, levels = c("Y", "X")))

# align the abundance dataset with the selected metadata
dat_kallisto_raw <- dat_kallisto_raw %>%
  select(rownames(metadata))

# prepare a limma/voom/edgeR/dream-compatible data structure
# edgeR doesn't handle mixed effects - dream does but the overall approach is different
dge <- DGEList(counts = dat_kallisto_raw,
               samples = metadata)
# calculate library size normalization factors
dge <- calcNormFactors(dge)

# define the model
form <- ~ study_arm + study_week + study_arm:study_week + (1|PID)

# estimate precision weights (mean-variance relationships) for each MAG. This approach builds on the 
# 'traditional' limma/voom approach to allow mixed effects. 
vobjDream <- voomWithDreamWeights(counts = dge, 
                                  formula = form, 
                                  data = metadata,
                                  span = 0.8, 
                                  plot = TRUE,
                                  save.plot = TRUE)

# Fit the models for each MAG, then perform hypothesis tests as defined in the contrast matrix (default or custom)
# degrees of freedom are estimated using the Satterthwaite approximation by default (Kenward-Roger is more accurate, but slower - recommended for < 10 samples)
fitmm <- dream(exprObj = vobjDream,
               formula = form, 
               data = metadata)

# exploring the modeling objects and contrasts - this can help you identify/specify your calls for results
L = getContrast(vobjDream, form, metadata, c('study_armX'))
plotContrasts(L)
# Examine design matrix
head(fitmm$design, 3)

# if you don't have a mixed model but still want to use dream, you'll need calls of the form below to perform the simpler testing
#fit_dream_mag_RUSF_baseline_quartile_model <- dream(exprObj = voom_dream_mag_RUSF_baseline_quartile_model, formula = form, data = anthro_rusf_baseline_quartiles, L)
#fit_dream_mag_RUSF_baseline_quartile_model <- eBayes(fit_dream_mag_RUSF_baseline_quartile_model)

# Extract results of DE testing for coefficients of interest
res.study_arm <- topTable(fitmm, coef='study_armX', adjust.method = "BH", n = Inf) %>%
  arrange(t)

res.interaction <- topTable(fitmm, coef='study_armX:study_week', adjust.method = "BH", n = Inf) %>%
  arrange(t)

# write results to file
write.table(res.study_arm, file = paste0(datestring, "_MDCF_POC_MAGs_dream_study_arm.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(res.interaction, file = paste0(datestring, "_MDCF_POC_MAGs_dream_interaction.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

# interpreting LME results using GSEA -------------------------------------

# each analysis will have a ranking strategy for the MAGs and a set of MAG groupings to be tested

# first, we'll look at MAGs ranked by their quartile association

# Gene sets are wlz-associated MAGs (negative and positive separately)

# load the WLZ association information
wlz_mags <- read.table("~/Box/MDCF_POC_Microbiome_Paper/MAGs/210614_0934_bin_1000_wlz_res_sig_vst.txt", header = TRUE, sep = "\t")

# create the sets
wlz_mag_sets <- split(wlz_mags$MAG, wlz_mags$sign)

# ranks are Q4 vs Q1 differences
ranks <- res.quartile$t
names(ranks) <- res.quartile$MAG

# run GSEA via fgsea
fgsea_res <- fgsea(wlz_mag_sets,
                   minSize = 10, # this may need tuning depending on your analysis
                   maxSize = length(ranks)/2, # this is a good starting point
                   stats = ranks,
                   eps = 0)

fgsea_res <- fgsea_res %>%
  arrange(desc(NES))
fgsea_res %>%
  filter(padj < 0.05)

# write results to file
fwrite(fgsea_res, file = paste0(datestring, "_MDCF_POC_MAGs_dream_quartile_fgsea_wlz_mags.txt"), sep = "\t", sep2 = c("", " ", ""))

# ranks are Q4 vs Q1 RATE differences
ranks <- res.interaction$t
names(ranks) <- res.interaction$MAG

fgsea_res <- fgsea(wlz_mag_sets,
                   minSize = 10,
                   maxSize = length(ranks)/2,
                   stats = ranks,
                   eps = 0)

fgsea_res <- fgsea_res %>%
  arrange(desc(NES))
fgsea_res %>%
  filter(padj < 0.05)

fwrite(fgsea_res, file = paste0(datestring, "_MDCF_POC_MAGs_dream_interaction_fgsea_wlz_mags.txt"), sep = "\t", sep2 = c("", " ", ""))

# checking directionality of enrichment results ---------------------------

# ...for "sanity"

# summarizing VST-normalized data
delta_POC <- dat_kallisto_vst.anthropometry %>%
  pivot_longer(c(-PID, -study_week, -SID, -gender, -enrollment_age, -baseline_age, -study_arm, -age_at_sampling, -study_phase_alt, -sample_number, -study_phase, -sample_type, -wlz), names_to = "bin", values_to = "ab") %>%
  pivot_wider(id_cols = c(study_arm, PID, bin), names_from = study_week, values_from = ab)

name.fix <- colnames(delta_POC)[!colnames(delta_POC) %in% c("study_arm", "PID", "bin")]
colnames(delta_POC)[!colnames(delta_POC) %in% c("study_arm", "PID", "bin")] <- paste0("wk_", name.fix)

# convert to long form
delta_POC.melt <- delta_POC %>%
  left_join(map.quartiles %>% select(-study_arm), by = "PID") %>%
  pivot_longer(c(-PID, -study_arm, -coef, -beta_quartile, -bin), names_to = "study_week", values_to = "ab")

# positive by quartile
dat.plot_pos <- delta_POC.melt %>%
  filter(bin == res.quartile$MAG[1],
         study_arm == "X",
         beta_quartile %in% c("Q1", "Q4")) %>%
  mutate(study_week = factor(study_week, levels = c("wk_0", "wk_2", "wk_4", "wk_8", "wk_12", "wk_16")))

ggplot(dat.plot_pos, aes(x = study_week, y = ab, group = beta_quartile, color = beta_quartile)) +
  geom_smooth(method = "lm")

# negative by arm
dat.plot_neg <- delta_POC.melt %>%
  filter(bin == res.quartile$MAG[nrow(res.quartile)],
         study_arm == "X",
         beta_quartile %in% c("Q1", "Q4")) %>%
  mutate(study_week = factor(study_week, levels = c("wk_0", "wk_2", "wk_4", "wk_8", "wk_12", "wk_16")))

ggplot(dat.plot_neg, aes(x = study_week, y = ab, group = beta_quartile, color = beta_quartile)) +
  geom_smooth(method = "lm")

# check the GSEA result - this function is from fgsea
plotEnrichment(wlz_mag_sets[["pos"]], ranks) + labs(title="pos")

plotEnrichment(wlz_mag_sets[["neg"]], ranks) + labs(title="neg")

# check the GSEA result - this function is from fgsea
plotGseaTable(wlz_mag_sets, 
              ranks, 
              fgsea_res, 
              gseaParam=0.5)

# fgsea by species

tax.filt <- tax %>%
  filter(MAG %in% lme_res_kallisto.full$MAG)

species_sets <- split(tax.filt$MAG, tax.filt$species)

# ranks are quartile associations
ranks <- limma_dream_MDCF2_quartile_model$t
names(ranks) <- limma_dream_MDCF2_quartile_model$MAG

fgsea_res <- fgsea(species_sets,
                   minSize = 5,
                   maxSize = length(ranks)/2,
                   stats = ranks,
                   eps = 0)

fgsea_res <- fgsea_res %>%
  arrange(desc(NES))
fgsea_res %>%
  filter(padj < 0.05)

tax.filt %>%
  filter(species == "s__Prevotella sp900290275") %>%
  select(MAG, genus, species)

sig.kallisto_vst %>%
  filter(species == "s__Prevotella sp900290275")
