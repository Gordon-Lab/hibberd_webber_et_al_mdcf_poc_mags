# rm(list=ls()) #Clear out old variables
setwd("C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/May2022_Fig4_follow_up_analysis")  # set working dir
library(fgsea)
library(tidyverse)
library(openxlsx)
library(ggrepel)
library(FactoMineR)
library(factoextra)
library(vegan) #adonis2

# Load beta-wlz quartiles:
quartile <- read.xlsx("C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/reanalysis_210826_onward/input/210810_1323_MDCF_POC_participant_beta_WLZ_quartiles.xlsx",
                      sheet=2)
colnames(quartile)
#[1] "PID"                     "study_arm"               "MAGabundance_model_coef" "delta_wlz"               "beta_quartile" 

#load vsd counts: Filters applied A) filtered out genes that were not part of 837 MAGs in Matt's abundance filtered set, 
# B) filtered out genes annotated as rRNA 
# C) set transcripts to zero for MAGs with abundance <.5 TPM 
# D) Filtered out genes with <10 raw counts total among all samples
  (load(file="210621_MDCF_POC_transcript_counts_vst_filtered_SumCtLT10_wTranscriptZeroing.Rda"))
  class(vsd)
#[1] "DESeqTransform"
#attr(,"package")
#[1] "DESeq2"
  vsd <- as.data.frame(assay(vsd))
  class(vsd)
  dim(vsd)
#[1] 1929056     350
View(head(vsd))
head(rownames(vsd))
#[1] "MBINA83BCE1B_00005" "MBINA83BCE1B_00010" "MBINA83BCE1B_00015" "MBINA83BCE1B_00020" "MBINA83BCE1B_00025" "MBINA83BCE1B_00030"
head(colnames(vsd))
#[1] "P01C5551102" "P01C5551105" "P01C5551107" "P01C6661102" "P01C6661105" "P01C6661107"

#load MAG abundance
MAG_vst <- read.table("C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/MAGs/210614_0934_MDCF_POC_MAG_1000_counts_vst_filt.txt", sep= "\t", header = TRUE, row.names = 1)

#Load a list of locus_tags with counts passing prevalence and abundance filter 
# Source: /lts/jglab/publications/MDCF_POC_Microbiome_Paper/RNAseq/210621_MDCF_POC_transcript_counts_filtered_SumCtLT10_wTranscriptZeroing_EdgeR_flt.tsv
# File description: kallito transcript raw counts filtered_wTranscriptZeroing and totoal counts <10, edgeR filtered at min.count = 5, min.prop = 0.70
raw_counts_flt <- read.table("C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/210621_MDCF_POC_transcript_counts_filtered_SumCtLT10_wTranscriptZeroing_EdgeR_flt.tsv", row.names=1, sep= "\t", header = TRUE)

fig4locus_tags <- c("MBINF3D6CCAE_00970", "MBINF3D6CCAE_07165", "MBINF3D6CCAE_01585","MBINCFF3C9CF_04015", "MBINF3D6CCAE_07470",
                    "MBINFE7566FD_09060", "MBINF3D6CCAE_11045", "MBINF3D6CCAE_06645", "MBIND2C09C28_07610", "MBIN0B7A9C00_12305", "MBIN7AA2F86E_01310")

#load the locus link info
(load("C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/Input/MAG_locus_link.Rda"))
#[1] "MAG_locus_link"
head(colnames(MAG_locus_link))
#[1] "MAG"   "locus"
head(MAG_locus_link$MAG)
#[1] "P01C555_maxbin.024.mod2"  "P01C555_maxbin.026.mod2"  "P01C555_maxbin.072.mod2" 
#[4] "P01C555_maxbin.077.mod2"  "P01C555_maxbin.081.mod2"  "P01C555_metabat.045.mod2"
head(MAG_locus_link$locus)
#[1] "MBINA83BCE1B" "MBIN31520C81" "MBIN0502F479" "MBIN8A183361" "MBINE5FAAB18" "MBIN66C42D2E"

#Load sample metadata
# load files "filesMeta"
(load("C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/reanalysis_210826_onward/input/sample_metadata.Rda"))
dim(filesMeta)
# 352  13
colnames(filesMeta)
# [1] "sample"          "PID"             "gender"          "enrollment_age"  "baseline_age"    "study_arm"       "age_at_sampling"
# [8] "study_week"      "sample_number"   "study_phase"     "study_phase_alt" "sample_type"     "path"       

#Left join quartile and metadata
filesMeta <- left_join(filesMeta, quartile)
map <- filesMeta
colnames(map)[1] <- "SID"

 #Load MAG data
MAG_meta_data <- read.xlsx("C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/reanalysis_210826_onward/input/210617_MDCF_POC_MAG_tables_rFriendly.xlsx", sheet=3)
  #Add on the MAG names
  vsd_tag <- separate(
  rownames_to_column(vsd, "locus_tag"),
  locus_tag,
  c("locus","tag"),
  sep = "_",
  remove = FALSE)
  
  vsd_tag <- left_join(vsd_tag, MAG_locus_link)
  View(head(vsd_tag))
  
#Make a locus link file
  locus_tag_MAG <- select(vsd_tag,c("locus_tag","locus","tag","MAG"))
  dim(locus_tag_MAG)
  #[1] 1929056       4
#Make a list of wlz-associated MAGs
  wlz_assoc_MAG_list <- MAG_meta_data[MAG_meta_data$MAG_abd_FDR_adjusted<.05,c("MAG")]
#Filter vst counts to include only WLZ associated MAGs (q<.05)
  vsd_wlz_assoc <- vsd_tag[vsd_tag$MAG %in% wlz_assoc_MAG_list,]
  dim(vsd_wlz_assoc)

#Reformat so that the df is numeric only
  rownames(vsd_wlz_assoc) <- vsd_wlz_assoc$locus_tag
  vsd_wlz_assoc <- select(vsd_wlz_assoc,!c("locus_tag","locus","tag","MAG"))

#Filter to match the final list of transcripts used in other analysis that were derived as follows:
  #######################################
  # Run PCA
  # run PCA via svd (missing values can be an issue here)
  # vst counts: Filters applied A) filtered out genes that were not part of 837 MAGs in Matt's abundance filtered set, 
  # B) filtered out genes annotated as rRNA
  # C) set transcripts to zero for MAGs with abundance <.5 TPM
  # D) Filtered out genes with <10 raw counts total among all samples
  # E) Filtered to include vst transcript counts for 222 WLZ associated MAGs
  # F) Filtered for edgeR with min.count = 5, min.prop = 0.70
  ################
  dim(raw_counts_flt)
  rm(vsd_wlz_assoc)
  vsd_wlz_assoc <- vsd[rownames(vsd) %in% rownames(raw_counts_flt),]
  dim(vsd_wlz_assoc)
  #[1] 27518   350
  dat_kallisto_vst <- vsd_wlz_assoc
  res_pca <- prcomp(t(dat_kallisto_vst), center = TRUE, scale = FALSE)

  # alternative PCA calculation tools
  #stats::princomp #spectral decomposition
  #stats::prcomp #svd
  #FactoMineR::PCA #svd
  
  # extract PC variances explained and make scree plot (both functions from factominer)
  pca_res.eig <- get_eigenvalue(res_pca)
  p <- fviz_eig(res_pca)
    pdf_filename <- "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/220808_transcript_pca_scree_plot.pdf"
    pdf(pdf_filename, width = 6, height = 6)
    print(p) #explicitly print the plot
    dev.off()
    
  # pca visualizations - these are messy and not terribly useful
  #fviz_pca_ind(res_pca)
  #fviz_pca_var(res_pca)
  
  # extract pc projections for samples - "ind" is individuals (samples), "var" is variables
  pca_res.ind <- get_pca_ind(res_pca)
  
  # add annotation metadata - this lets me plot outside factominer with more customization
  pca_res.ind.anno <- pca_res.ind$coord %>%
    as.data.frame() %>%
    rownames_to_column("SID") %>%
    left_join(map %>% select(SID, PID, study_week, study_arm, delta_wlz, beta_quartile, study_phase_alt), by = "SID")

  # adonis (permanova)
  ########################
  pca_res.ind.sub <- pca_res.ind.anno %>%
    #  filter(study_week == "0") %>%
    select(SID, Dim.1, Dim.2) %>%
    column_to_rownames("SID")

  # permutations not blocked
  adonis_res <- adonis2(pca_res.ind.sub ~ study_arm*study_week, data = pca_res.ind.anno, method = "euclidean")
  adonis_res
  write.table(adonis_res, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/220727_transcript_pca_PERMANOVA_not_blocked.txt", sep = "\t", row.names = FALSE)
  
  # permutations are blocked within participants, but study_arm doesn't vary within participants, so I'm not sure this is appropriate
  perm <- how(nperm = 199)
  setBlocks(perm) <- with(pca_res.ind.anno, PID)
  adonis_res_blocked <- adonis2(pca_res.ind.sub ~ study_arm*study_week, data = pca_res.ind.anno, method = "euclidean", permutations = perm)
  adonis_res_blocked
  write.table(adonis_res_blocked, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/220727_transcript_pca_PERMANOVA_blocked.txt", sep = "\t", row.names = FALSE)
  
  
  # checking pairwise by week
  pca_res.ind.anno.sub <- pca_res.ind.anno %>%
    filter(study_week == 12) %>%
    select(SID, Dim.1, Dim.2, study_arm, study_week, study_phase_alt, beta_quartile) %>%
    column_to_rownames("SID")
  
  # permutations not blocked
  adonis_res_wk12 <- adonis2(pca_res.ind.anno.sub %>% select(Dim.1, Dim.2) ~ study_arm, data = pca_res.ind.anno.sub, method = "euclidean")
  adonis_res_wk12
  write.table(adonis_res_wk12, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/220727_transcript_pca_PERMANOVA_wk12_not_blocked.txt", sep = "\t", row.names = FALSE)
  
  # calculating centroids and plotting with ggplot - more flexible
  
  pca_res.ind.anno.cent <- pca_res.ind.anno %>%
    #  filter(study_arm != "RUSF") %>%
    group_by(study_arm, study_phase_alt, study_week) %>%
    summarize_at(vars(starts_with("Dim")), mean)
  
  # The palette with black:
  pca_res.ind.anno$study_week <- as.factor(pca_res.ind.anno$study_week)
  pca_res.ind.anno.cent$study_week <- as.factor(pca_res.ind.anno.cent$study_week)
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  p <- ggplot(pca_res.ind.anno, aes(x = Dim.1, y = Dim.2, color = study_phase_alt, shape = study_arm)) +
    geom_point(size = 2) +
    scale_color_manual(values = cbbPalette) +
    coord_equal() +
    theme_bw() +
    xlab(paste0(rownames(pca_res.eig[1,]), " ", round(pca_res.eig[1,]$variance.percent, 2), "% variance")) +
    ylab(paste0(rownames(pca_res.eig[2,]), " ", round(pca_res.eig[2,]$variance.percent, 2), "% variance"))
  pdf_filename <- "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/220808_transcript_pca_all_samples.pdf"
  pdf(pdf_filename, width = 8, height = 6)
  print(p) #explicitly print the plot
  dev.off()
  
  p <- ggplot(pca_res.ind.anno %>% filter(study_phase_alt %in% c("baseline", "treatment_phase3")), aes(x = Dim.1, y = Dim.2, group = paste(study_week, study_arm), color = study_phase_alt, shape = study_arm)) +
    geom_point(size = 2) +
    scale_color_manual(values = cbbPalette) +
    stat_ellipse() +
    theme_bw() +
    coord_equal() +
    xlab(paste0(rownames(pca_res.eig[1,]), " ", round(pca_res.eig[1,]$variance.percent, 2), "% variance")) +
    ylab(paste0(rownames(pca_res.eig[2,]), " ", round(pca_res.eig[2,]$variance.percent, 2), "% variance"))
  pdf_filename <- "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/220808_transcript_pca_baseline-endpoint_wVariance.pdf"
  pdf(pdf_filename, width = 8, height = 6)
  print(p) #explicitly print the plot
  dev.off()
  
  p <- ggplot(pca_res.ind.anno.cent, aes(x = Dim.1, y = Dim.2, color = study_phase_alt, shape = study_arm)) +
    geom_point(size = 4) +
    geom_path(aes(group = study_arm)) +
    scale_color_manual(values = cbbPalette) +
    coord_equal() +
    theme_bw() +
    xlab(paste0(rownames(pca_res.eig[1,]), " ", round(pca_res.eig[1,]$variance.percent, 2), "% variance")) +
    ylab(paste0(rownames(pca_res.eig[2,]), " ", round(pca_res.eig[2,]$variance.percent, 2), "% variance"))
  #  guides(color = FALSE)
  pdf_filename <- "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/220808_transcript_pca_centroids_time_arm.pdf"
  pdf(pdf_filename, width = 8, height = 6)
  print(p) #explicitly print the plot
  dev.off()
  
  pca_res.ind.anno.cent2 <- pca_res.ind.anno %>%
    #  filter(study_arm != "Y") %>%
    group_by(study_arm, beta_quartile, study_phase_alt) %>%
    summarize_at(vars(starts_with("Dim")), mean)
  
  p <- ggplot(pca_res.ind.anno.cent2, aes(x = Dim.1, y = Dim.2, color = study_phase_alt, shape = study_arm)) +
    geom_point(size = 4) +
    geom_path(aes(group = study_arm)) +
    scale_color_manual(values = cbbPalette) +
    theme_bw() +
    coord_equal() +
    xlab(paste0(rownames(pca_res.eig[1,]), " ", round(pca_res.eig[1,]$variance.percent, 2), "% variance")) +
    ylab(paste0(rownames(pca_res.eig[2,]), " ", round(pca_res.eig[2,]$variance.percent, 2), "% variance")) +
    facet_grid(study_arm ~ beta_quartile) 
  pdf_filename <- "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/220808_transcript_pca_centroids_time_arm_quartile.pdf"
  pdf(pdf_filename, width = 12, height = 4)
  print(p) #explicitly print the plot
  dev.off()
  
  p <- ggplot(pca_res.ind.anno.sub, aes(x = Dim.1, y = Dim.2, color = beta_quartile, shape = study_arm)) +
    geom_point(size = 4) +
    scale_color_manual(values = cbbPalette) +
    stat_ellipse() +
    theme_bw() +
    coord_equal() +
    xlab(paste0(rownames(pca_res.eig[1,]), " ", round(pca_res.eig[1,]$variance.percent, 2), "% variance")) +
    ylab(paste0(rownames(pca_res.eig[2,]), " ", round(pca_res.eig[2,]$variance.percent, 2), "% variance"))
  pdf_filename <- "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/220808_transcript_pca_samples_wVariance_arm_quartile.pdf"
  pdf(pdf_filename, width = 8, height = 6)
  print(p) #explicitly print the plot
  dev.off()
  
  # permutations not blocked
  adonis_res_quartile_both_arms <- adonis2(pca_res.ind.anno.sub %>% select(Dim.1, Dim.2) ~ beta_quartile, data = pca_res.ind.anno.sub, method = "euclidean")
  adonis_res_quartile_both_arms
  write.table(adonis_res_quartile_both_arms, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/220727_transcript_pca_PERMANOVA_quartile_both_arms_not_blocked.txt", sep = "\t", row.names = FALSE)
  
  
  # mag (variable) contributions to coords
  
  #mag_number_key <- read_excel("220725 MDCF POC PM mag numbering key.xlsx")
  
  #wlz_mag_ids <- sig.kallisto_vst %>%
  #  inner_join(mag_number_key, by = "MAG") %>%
  #  select(MAG_ID, sign)
  
  #target_mags <- read_excel("~/Documents/Projects/human_studies/MDCF POC/220422 MAG abundance by quartile re-analysis/220422 Fig4 target mags.xlsx")
  
  tax <- MAG_meta_data %>% dplyr::rename(genus = Genus) %>% dplyr::rename(species = Species)
  tax$sign <- sign(tax$beta)
  pca_res.var <- get_pca_var(res_pca)
  pca_res.var.anno <- pca_res.var$contrib %>%
    as.data.frame() %>%
    rownames_to_column("locus_tag") %>%
    left_join(locus_tag_MAG, by = "locus_tag") %>%
    left_join(tax %>% select(MAG, genus, species, MAG_ID), by = "MAG") %>%
    mutate(label = paste(locus_tag, species, sep = " | "),
           label = factor(label, levels = label))
  
  write.table(pca_res.var.anno %>% select(locus_tag, MAG_ID, MAG, genus, species, Dim.1, Dim.2, Dim.3), "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/220727_transcript_pca_contributions.txt", sep = "\t", row.names = FALSE)
  
  # contributions to Dim 1
  Dim.1 <- pca_res.var.anno %>% 
    select(MAG_ID, Dim.1, species) %>%
    slice_max(Dim.1, n = 100)
  Dim.1
  #write.table(Dim.1, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/220727_PC1_max_contributors.txt", sep = "\t", row.names = FALSE)
  
  
  pca_res.var.anno %>% 
    select(MAG_ID, Dim.1, genus, species) %>%
    arrange(desc(Dim.1)) %>%
    mutate(rank = row_number()) %>%
    filter(genus == "Prevotella")
  
  pca_res.var.anno %>% 
    select(MAG_ID, Dim.1, species) %>%
    arrange(desc(Dim.1)) %>%
    mutate(rank = row_number()) %>%
    filter(MAG_ID %in% c("MAG_0018", "MAG_0019"))
  
  pca_res.var.anno %>% 
    select(locus_tag, MAG_ID, Dim.1, species) %>%
    arrange(desc(Dim.1)) %>%
    mutate(rank = row_number()) %>%
    filter(locus_tag %in% fig4locus_tags)
  
  # contributions to Dim 2
  pca_res.var.anno %>% 
    select(MAG_ID, Dim.2, species) %>%
    slice_max(Dim.2, n = 10)
  
  pca_res.var.anno %>% 
    select(MAG_ID, Dim.2, genus, species) %>%
    arrange(desc(Dim.2)) %>%
    mutate(rank = row_number()) %>%
    filter(genus == "Bifidobacterium")
  
  pca_res.var.anno %>% 
    select(MAG_ID, Dim.2, species) %>%
    arrange(desc(Dim.2)) %>%
    mutate(rank = row_number()) %>%
    filter(MAG_ID %in% c("MAG_0018", "MAG_0019"))
  
  pca_res.var.anno %>% 
    select(locus_tag, MAG_ID, Dim.2, species) %>%
    arrange(desc(Dim.2)) %>%
    mutate(rank = row_number()) %>%
    filter(locus_tag %in% fig4locus_tags)
  
  # matrix shuffling - this helps you decide how many PCs capture meaningful quantities of variation
  
  elements <- ncol(dat_kallisto_vst)*nrow(dat_kallisto_vst)
  dat_kallisto_vst.rand <- matrix(sample(c(t(dat_kallisto_vst)), elements, replace = FALSE), ncol = ncol(dat_kallisto_vst)) %>%
    as.data.frame()
  colnames(dat_kallisto_vst.rand) <- colnames(dat_kallisto_vst)
  rownames(dat_kallisto_vst.rand) <- rownames(dat_kallisto_vst)
  res_pca.rand <- prcomp(t(dat_kallisto_vst.rand), center = TRUE, scale = FALSE)
  noise_threshold <- get_eigenvalue(res_pca.rand)[1,]$variance.percent
  
  n_meaningful_PCs <- pca_res.eig %>%
    filter(variance.percent > noise_threshold) %>% nrow()
  meaningful_PCs <- pca_res.eig[1:meaningful_PCs,]
  meaningful_PCs
  
  # checking directionality of variables along PCs
  
  pca_res.dir <- pca_res.var$coord %>%
    as.data.frame() %>%
    rownames_to_column("locus_tag") %>%
    left_join(locus_tag_MAG, by = "locus_tag") %>%
    select(locus_tag, MAG, Dim.1, Dim.2, Dim.3) %>%
    left_join(tax %>% select(MAG, MAG_ID, genus, species), by = "MAG")
   
  write.table(pca_res.dir %>% select(locus_tag, MAG_ID, MAG, Dim.1, Dim.2, Dim.3, species), "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/220727_transcript_pca_directionality.txt", sep = "\t", row.names = FALSE)
  
  pca_res.dir %>%
    slice_max(Dim.1, n = 10)
  pca_res.dir %>%
    slice_min(Dim.1, n = 10)
  
  pca_res.dir %>%
    slice_max(Dim.2, n = 10)
  pca_res.dir %>%
    slice_min(Dim.2, n = 10)
  
  # running fgsea to check for species or wlz enrichment along a pc
  fgsea_res.bulk9
  
  sets <- split(pca_res.dir$locus_tag, pca_res.dir$species)
  sets <- split(pca_res.dir$locus_tag, pca_res.dir$genus)
  sets <- split(tax$locus_tag, tax$sign)
  
  ranks <- pca_res.dir$Dim.1
  names(ranks) <- pca_res.dir$locus_tag
  ranks <- sort(ranks)
  
  fgsea_res <- fgsea(pathways = sets,
                     stats = ranks,
                     eps = 0) %>%
    filter(padj < 0.05) %>%
    arrange(desc(NES))
  View(fgsea_res)
  
    fgsea_res.bulk <- fgsea_res.bulk %>%
    bind_rows(fgsea_res %>%
                mutate(dim = "Dim.1") %>%
                rowwise() %>%
                mutate(leadingEdge = paste(unlist(leadingEdge), collapse = ",")) %>%
                select(dim, set = pathway, padj, NES, size, leadingEdge))
  
  write.table(fgsea_res.bulk, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/220728_species_PCA_wlz_species_enrichment.txt", sep = "\t", row.names = FALSE)
  