#EdgeR info: https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html
#ptmixed: https://cran.r-project.org/web/packages/ptmixed/vignettes/Overview_functionalities_ptmixed.html
#2020 Ptmixed paper: https://arxiv.org/pdf/2004.11193.pdf
# Info on NBmixed: https://github.com/rtsonaka/NBmixed_RNAseq
# Dream vignette: https://bioconductor.org/packages/release/bioc/vignettes/variancePartition/inst/doc/dream.html
# Limma/vooom vignett: https://f1000research.com/articles/5-1408/v3
#Clear out old variables
rm(list=ls()) 
# library(stringr) # str_sub function for getting substrings

# library(ggpubr)
# library(broom)
library(openxlsx)
library(ggplot2)
library(ggsignif)
library(pheatmap) #For heatmaps
library(tidyverse)
library('variancePartition')
library('tximport')
library('BiocParallel')
library('edgeR')
library('DESeq2')
library('RColorBrewer')
library(superSeq) #For subsampeling 
library(subSeq)   #For subsampeling
library(Biobase)
library("vsn")
library("Glimma")
library("fgsea")
library(data.table) # Needed for GSEA output (fwrite)
library(ggpubr)
library(janitor)
# library("ptmixed")
# source("/scratch/jglab/dmwebber/R/NBmixed_RNAseq-master/NBmixed_RNAseq-master/MainFunction.R")       
# source("/scratch/jglab/dmwebber/R/NBmixed_RNAseq-master/NBmixed_RNAseq-master/Simulate_Data.R")       
#Source above script for NBmixed: https://github.com/rtsonaka/NBmixed_RNAseq 

## 1. User input
##############

setwd("/scratch/jglab/dmwebber/RNAseq/R_analysis")

# Read MAG taxa and wlz association
MAG_taxa <- read.xlsx("/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/input/210617_MDCF_POC_MAG_tables_rFriendly.xlsx", sheet = 3, startRow = 1, colNames = T, rowNames = F)

##################################################### 
#Load MAG-locus_tag link                            #
load(file="/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/input/MAG_locus_link.Rda")      # object: MAG_locus_link
#####################################################
MAG_taxa <- left_join(MAG_locus_link, MAG_taxa)
  # Change the MAG_ID to reflect current convention
  MAG_taxa <- MAG_taxa %>% 
    mutate(MAG_ID = str_replace(MAG_ID, "MAG_", "MAGbg"))
  MAG_taxa$MAG_ID_species <- paste(MAG_taxa$MAG_ID,MAG_taxa$Species,sep="_")

###### Load counts
load(file="/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/input/210621_MDCF_POC_transcript_counts_filtered_SumCtLT10_wTranscriptZeroing_EdgeR_flt.Rda")    # object: counts_EdgeR_flt

##### Load meta data
load(file="/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/input/sample_metadata.Rda")  # object: filesMeta

# Input quartiles
quartiles <- read.xlsx("/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/input/210810_1323_MDCF_POC_participant_beta_WLZ_quartiles.xlsx", sheet = 2, startRow = 1, colNames = T, rowNames = F)
# Merge quartiles with metatdata
filesMeta <- left_join(filesMeta, quartiles)
t12_SID <- filesMeta[filesMeta$study_week=="12" & filesMeta$study_arm =="X",c("sample")] # limit the vsd file to t12 and MDCF2


# load(file="./1000MAGs/out/annotations_062621.Rda")  # object: annotations
TableS5 <- read.xlsx("/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/input/TableS5_210819_r.xlsx", sheet = 1, startRow = 1, colNames = T, rowNames = F)
TableS5_all_caz <- TableS5[!is.na(TableS5$CAZyme),] #Filter to all CAZyme annotations
dim(TableS5_all_caz)
TableS5_mcseed_only <- TableS5[!TableS5$Involved=="CAZyme",]
dim(TableS5_mcseed_only)
# Split the Phenotype column by ;
TableS5_mcseed_only <- separate_rows(TableS5_mcseed_only, Phenotype, sep = "; ",  convert = TRUE)  
TableS6 <- read.xlsx("/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/input/TableS6_r.xlsx", sheet = 1, startRow = 1, colNames = T, rowNames = F)
glossary <- read.xlsx("/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/input/TableS6_r.xlsx", sheet = 2, startRow = 1, colNames = T, rowNames = F)
# Convert TableS6 (wide) to tall format
TableS6_tall <- pivot_longer(TableS6, cols=2:ncol(TableS6), names_to = "Phenotype", values_to ="Phenotype_pres")
dim(TableS6_tall)
TableS5_mcseed_only <- left_join(TableS5_mcseed_only, TableS6_tall)
dim(TableS5_mcseed_only)
# Check for NA is the merged "Phenotype_pres" col
sum(is.na(TableS5_mcseed_only$Phenotype_pres))

#Add MAG_ID to TableS5
TableS5_mcseed_only <- left_join(TableS5_mcseed_only, MAG_taxa[,c("MAG","MAG_ID", "Genus", "Species")])

# Output the new S5, which includes a columns for "Phenotype_pres" to indicate wither the binary penotype is 1 (present) or 0 (absent) for an annotated gene. 
write.csv(TableS5_mcseed_only, "/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/TableS5mod_210830.csv")

#Filter S5 based on presence (1) or absence (0) of predicted phenotype in MAG with positive (1) binary phenotype in S6
TableS5_mcseed_only_BP1 <- TableS5_mcseed_only[TableS5_mcseed_only$Phenotype_pres==1,]
#Annotated genes with MAG BP=1 
nrow(TableS5_mcseed_only_BP1)
#Genes filtered out d/t MAM BP=0
sum(TableS5_mcseed_only$Phenotype_pres==0)


# Load glycan data (objects: "linkage_delta_wk12" "mono_delta_wk12" "vst_delta_wk12")
load(file="/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/input/210802_fecal_carb_vs_vst_transcripts_for_cyrus.RData")
#Fix duplicate colname (replace "2_X_X_Hexose_II" with "2_X_X_Hexose_I")
colnames(linkage_delta_wk12)[length(colnames(linkage_delta_wk12))-1] <- "2_X_X_Hexose_I"
colnames(linkage_delta_wk12) <- janitor::make_clean_names(colnames(linkage_delta_wk12),case = "none")

# Load vst transcripts
load(file= "/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/input/210621_MDCF_POC_transcript_counts_vst_filtered_SumCtLT10_wTranscriptZeroing.Rda") #Object: vsd
vst_counts <- as.data.frame(assay(vsd))
vst_counts <- t(vst_counts)
vst_counts_wk12 <- vst_counts[rownames(vst_counts) %in% t12_SID, colnames(vst_counts) %in% colnames(vst_delta_wk12)] #Limt the vsd counts dataset to those transcripts in the vsd df
rownames(vst_counts_wk12) <- strtrim(rownames(vst_counts_wk12),7)

#List of monos to exclude
monosac_above_LOD <- c("Arabinose", "Fructose",  "Fucose",  "GalA", "Galactose","Glucose","Mannose", "Rhamnose","Ribose","Xylose")
# Read in GLYCAN data:
fec_monosac <- read.xlsx("/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/input/2021_07_02_Monosaccharide and Linkage analysis of POC study and Diets_20210728.xlsx", sheet = 1, startRow = 2, colNames = T, rowNames = F)
rownames(fec_monosac) <- fec_monosac$Sample.ID
# Make a list of the duplicated IDs  
fec_monosac_rep_ID <- c("P13C5551107_01", "P13C5551107_02", "P13C5551107_03", "P25C7771105_01", "P25C7771105_02", "P57C7771107_01", "P57C7771107_02", "P57C7771107_03")
# Exclude replicate samples
fec_monosac <- fec_monosac[!fec_monosac$Sample.ID %in% fec_monosac_rep_ID,]
fec_monosac$Sample.number <- NULL
fec_monosac <- fec_monosac[,colnames(fec_monosac) %in% monosac_above_LOD]

fec_linkage <- read.xlsx("/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/input/2021_07_02_Monosaccharide and Linkage analysis of POC study and Diets_20210728.xlsx", sheet = 5, startRow = 2, colNames = T, rowNames = F)
colnames(fec_linkage) <- janitor::make_clean_names(colnames(fec_linkage),case = "none")
rownames(fec_linkage) <- fec_linkage$Sample_ID
fec_linkage$Sample_ID <- NULL
fec_linkage$Sample_number <- NULL
# Make a list of which IDs are replicates  
fec_linkage_rep_ID <- c("P13C5551107_01", "P13C5551107_02", "P13C5551107_03", "P25C7771105_01", "P25C7771105_02", "P57C7771107_01", "P57C7771107_02", "P57C7771107_03")
# Exclude replicate samples
fec_linkage <- fec_linkage[!rownames(fec_linkage) %in% fec_linkage_rep_ID,]

#Check to see if the linkage names for Rob's set are in the same format as those imported
linkage_delta_wk12[!colnames(linkage_delta_wk12) %in% colnames(fec_linkage)] #Display any columns from Rob's set that are not part of my set
colnames(fec_linkage)[!colnames(fec_linkage) %in% colnames(linkage_delta_wk12)] #Display any columns from my set that are not part of Rob's set

# Read in the food glycan data:
# Read in GLYCAN data:
food_monosac <- read.xlsx("/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/input/2021_07_02_Monosaccharide and Linkage analysis of POC study and Diets_20210728.xlsx", sheet = 3, startRow = 1, colNames = T, rowNames = T)
food_linkage <- read.xlsx("/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/input/2021_07_02_Monosaccharide and Linkage analysis of POC study and Diets_20210728.xlsx", sheet = 7, startRow = 1, colNames = T, rowNames = T)
colnames(food_linkage) <- janitor::make_clean_names(colnames(food_linkage),case = "none")
#Check to see if the linkage names for Rob's set are in the same format as those imported
colnames(food_linkage)[!colnames(food_linkage) %in% colnames(fec_linkage)] #Display any columns from the food glycan set that are not part of the fecal glycans

# Load food freq questionnaire data 
FFQ <- read.csv(file = '/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/FFQ.csv', header = TRUE)
#################################################
#################################################
#################################################
# Make sure that the metadata matches with the counts
filesMeta <- filesMeta[filesMeta$sample %in% colnames(counts_EdgeR_flt),]
all.equal(filesMeta$sample, colnames(counts_EdgeR_flt))

#Make study_arm a factor
filesMeta$study_arm <- as.factor(filesMeta$study_arm)

#Make sure that the control is the reference level:
filesMeta$study_arm <- relevel( filesMeta$study_arm, ref="Y")
filesMeta$study_arm

#Make a DGEList in EdgeR    #optional: counts[keep.exprs,]
x = DGEList(counts=counts_EdgeR_flt, group=filesMeta$study_arm)

#Normalize the data again since library size has changed:
x <- calcNormFactors(x, method = "TMM") #method = "TMMwsp", "TMM"

# Center/scale study week
filesMeta$study_week_nonSl <- filesMeta$study_week 
#filesMeta$study_week <- log2(filesMeta$study_week + abs(min(filesMeta$study_week))+1)
filesMeta$study_week <- scale(filesMeta$study_week, center = TRUE, scale = FALSE)

######################################### 
#Model C: transcript ~ study_arm + study_week study_arm*study_week
#####################
# Setup model matrix                            
design.mat2 <- model.matrix(~ filesMeta$study_arm + filesMeta$study_week + filesMeta$study_arm:filesMeta$study_week)
colnames(design.mat2) <- gsub("filesMeta\\$", "", colnames(design.mat2)) #Reformat the colnames using gsub. Don't forget to escape $ with double \\
design.mat2
# A. Fit the common dispersion with edgeR
d2c <- estimateGLMCommonDisp(x,design.mat2)
# B. Fit the tagwise dispersion with edgeR
d2c <- estimateGLMTrendedDisp(d2c, design.mat2, method="power")
# C. Estimate a generalized linear model (glm) fit with edgeR
d2c <- estimateGLMTagwiseDisp(d2c,design.mat2)
pdf("/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/arm_week_ArmWeek/EdgeR_BCV_GLMTrendedDisp_contVars.pdf", width = 7, height = 7)
plotBCV(d2c)
dev.off()
fit2 <- glmQLFit(d2c, design.mat2)

# View fit:
View(fit2)
# View fit coeffiecients
View(fit2$coefficients)
# study_armX
coef2_fit2 <- glmQLFTest(fit2, coef=2)
study_armX_topTags <- topTags(coef2_fit2, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = .05) #Coef: study_armX
study_armX_moreTags <- topTags(coef2_fit2, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = 1) #Coef: study_armX
write.csv(study_armX_topTags, "/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/arm_week_ArmWeek/study_armX_topTags.csv") 
sig_study_armX_topTags <- decideTestsDGE(coef2_fit2, adjust.method="BH", p.value=0.05)
summary(sig_study_armX_topTags)  #517 down, 788 up (344d 2726u)

# study_week:
coef3_fit2 <- glmQLFTest(fit2, coef=3)
study_week_topTags <- topTags(coef3_fit2, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = .05) #Coef: study_week
study_week_moreTags <- topTags(coef3_fit2, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = 1) #Coef: study_week
write.csv(study_week_topTags, "/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/arm_week_ArmWeek/study_week_topTags.csv") # 7 DE genes
sig_study_week_topTags <- decideTestsDGE(coef3_fit2, adjust.method="BH", p.value=0.05)
summary(sig_study_week_topTags) # 66 down, 232 up(118d, 560u)


#study_armX:study_week
coef4_fit2 <- glmQLFTest(fit2, coef=4)
study_arm_week_topTags <- topTags(coef4_fit2, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = .05) #Coef: study_armX:week
study_arm_week_moreTags <- topTags(coef4_fit2, n = Inf, adjust.method = "BH", sort.by = "PValue", p.value = 1) #Coef: study_armX:week
write.csv(study_arm_week_topTags, "/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/arm_week_ArmWeek/1000MAGs/out/wlz_topTags.csv") # 64 DE genes
sig_study_arm_week_topTags <- decideTestsDGE(coef4_fit2, adjust.method="BH", p.value=0.05)
summary(sig_study_arm_week_topTags) # 136 down, 37 up (271d, 171u)

###############  
##########
## GSEA ##
##########
###############

# Prep DG expression data for GSEA
################
# Model C   #Model C: transcript ~ study_arm + study_week + study_arm*study_week
################
study_armX_moreTags_df <- rownames_to_column(as.data.frame(study_armX_moreTags), var="locus_tag")
study_week_moreTags_df <- rownames_to_column(as.data.frame(study_week_moreTags), var="locus_tag")
study_arm_week_moreTags_df <- rownames_to_column(as.data.frame(study_arm_week_moreTags), var="locus_tag")
ModelC_moreTags_df <- bind_rows(list(study_armX = study_armX_moreTags_df, study_week = study_week_moreTags_df, arm.week = study_arm_week_moreTags_df), .id = "group")
dim(ModelC_moreTags_df[!duplicated(ModelC_moreTags_df$locus_tag)  ,]) #Check to see that we still have the apropriate number of lines
ModelC_TopTags_df <- ModelC_moreTags_df[ModelC_moreTags_df$FDR<.1, ] # Get top significant tags (27531)
save(ModelC_moreTags_df, file = "/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/arm_week_ArmWeek/all_tags_allMAGs_EdgeR_model_Arm_Week_ArmbyWeek.Rda") #Save object: ModelC_moreTags_df
write.csv(ModelC_moreTags_df, "/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/arm_week_ArmWeek/all_tags_allMAGs_EdgeR_model_Arm_Week_ArmbyWeek.csv")
#########
# Load: ModelC_moreTags_df
load(file= "/scratch/jglab/dmwebber/RNAseq/R_analysis/1000MAGs/out/WlzMAGs_ContVars/model_Arm_Week_ArmbyWeek/all_tags_allMAGs_EdgeR_model_Arm_Week_ArmbyWeek_210710.Rda") #Load object: ModelC_moreTags_df
#########

# Merge tags with annotations then filter down to wlz-associated MAGs
# ModelC_moreTags_df
ModelC_MoreTags_df_wlz <- separate(
  ModelC_moreTags_df,
  locus_tag,
  c("locus","tag"),
  sep = "_",
  remove = FALSE)
ModelC_MoreTags_df_wlz <- left_join(ModelC_MoreTags_df_wlz, MAG_locus_link)
ModelC_MoreTags_df_wlz <- left_join(ModelC_MoreTags_df_wlz, MAG_taxa[, c("MAG","MAG_abd_FDR_adjusted")])
dim(ModelC_MoreTags_df_wlz[!duplicated(ModelC_MoreTags_df_wlz$locus_tag),]) #Make sure that we still have only 27531 unique locus tags
ModelC_MoreTags_df_wlz <- ModelC_MoreTags_df_wlz[ModelC_MoreTags_df_wlz$MAG_abd_FDR_adjusted<0.05,] #Keep only the transcripts from the wlz-associated MAGs
ModelC_MoreTags_df_wlz$MAG_abd_FDR_adjusted <- NULL 
dim(ModelC_MoreTags_df_wlz[!duplicated(ModelC_MoreTags_df_wlz$locus_tag),]) #Make sure that we now have 7747 unique locus tags that are from wlz-associated MAGs
save(ModelC_MoreTags_df_wlz, file = "/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/arm_week_ArmWeek/all_tags_wlaMAGs_EdgeR_model_Arm_Week_ArmbyWeek.Rda") #Save object: ModelC_moreTags_df
write.csv(ModelC_MoreTags_df_wlz, "/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/arm_week_ArmWeek/all_tags_wlzMAGs_EdgeR_model_Arm_Week_ArmbyWeek.csv")
my.contrasts <- c("study_armX","study_week","arm.week")

more_tags_EdgeR <- ModelC_MoreTags_df_wlz
# Make -log10(PValue)
more_tags_EdgeR$log10PValue <- -log10(more_tags_EdgeR$PValue)
more_tags_EdgeR[more_tags_EdgeR$logFC<0, c("log10PValue")] <- -(abs(more_tags_EdgeR[more_tags_EdgeR$logFC<0, c("log10PValue")]))  #When logFC<0 convert log10PValue to negative (abs function added for protection)
######################################
seq <- seq(from=1, to=length(my.contrasts), by=1)
rm(gsea_df,df)
gsea_df= data.frame()
for(i in seq){
  #Clear old vars
  rm(more_tags_EdgeR_ranks, df, more_tags_EdgeR_i)
  contrast_name <- (my.contrasts)[i]
  # Filter to contrast of interst(e.g.X12vsY12)
  more_tags_EdgeR_i <- more_tags_EdgeR[more_tags_EdgeR$group==contrast_name,]
  #Apply the rownames
  rownames(more_tags_EdgeR_i) <- more_tags_EdgeR_i$locus_tag
  # Sort decending by FDR
  more_tags_EdgeR_i <- more_tags_EdgeR_i[order(more_tags_EdgeR_i$log10PValue, decreasing = TRUE),] #Sort by log10PValue
  #Get ranks
  more_tags_EdgeR_ranks <- more_tags_EdgeR_i[,c("log10PValue")]
  #apply rowlabels
  names(more_tags_EdgeR_ranks) <- more_tags_EdgeR_i$locus_tag
  # filter the annotation list to those transcripts that are being tested
  # Use new phenotype annotions: "TableS5_mcseed_only_BP1"
  
  annotations_keep <- TableS5_mcseed_only_BP1[TableS5_mcseed_only_BP1$locus_tag %in% more_tags_EdgeR$locus_tag,]
  annotations_keep$MAG_species <- paste(annotations_keep$MAG_ID, annotations_keep$Species, sep = "_" )
  
  # Setup pathways
  #mcSEED "Phenotypes"
  phenotype_keep <- split(annotations_keep$locus_tag, as.factor(annotations_keep$Phenotype))
  
  #MAGs
  MAGPaths_keep <- split(annotations_keep$locus_tag, as.factor(annotations_keep$MAG_species)) 
  #Genus
  GenusPaths_keep <- split(annotations_keep$locus_tag, as.factor(annotations_keep$Genus))
  ################
  
  #GSEA 
  #mcSEED "Phenotypes"
  fgsea_phenotype_kp <- fgsea(phenotype_keep, more_tags_EdgeR_ranks, minSize = 10, eps=0)
  fgsea_phenotype_kp <- cbind(contrast = contrast_name, path = "mcSEED_phenotype", fgsea_phenotype_kp)  
  
  # MAGs
  fgsea_MAGs_kp <- fgsea(MAGPaths_keep, more_tags_EdgeR_ranks, minSize = 10, eps=0)
  fgsea_MAGs_kp <- cbind(contrast = contrast_name, path = "MAGs", fgsea_MAGs_kp)
  
  # Genus
  fgsea_Genus_kp <- fgsea(GenusPaths_keep, more_tags_EdgeR_ranks, minSize = 10, eps=0)  #nPermSimple = 10000
  fgsea_Genus_kp <- cbind(contrast = contrast_name, path = "Genus", fgsea_Genus_kp) 
  
  #concatenate results of GSEA
  df <- rbind(fgsea_phenotype_kp, fgsea_MAGs_kp, fgsea_Genus_kp)
  df_rows <- nrow(df)
  
  if(is.null(gsea_df)){ gsea_df_rows = 0} else {gsea_df_rows <- nrow(gsea_df)}
  gsea_df[(gsea_df_rows+1):(gsea_df_rows+df_rows),1:10] <- df
}

############## 
# Model C output: transcript ~ study_arm + study_week + study_arm*study_week
# Output GSEA results with LeadingEdge
############
fwrite(gsea_df, file="/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/arm_week_ArmWeek/GSEA_wlzMAGs_all_results.txt", sep="\t", sep2=c("", " ", ""))
save(gsea_df, file="/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/arm_week_ArmWeek/GSEA_wlzMAGs_all_results.Rda")

# Filter to significant results (Padj < 0.1)
gsea_df_sig <- gsea_df[gsea_df$padj<0.1 & !is.na(gsea_df$padj), ]
gsea_df_sig <- left_join(gsea_df_sig, glossary, by = c("pathway" = "Phenotype"))
fwrite(gsea_df_sig, file="/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/arm_week_ArmWeek/GSEA_wlzMAGs_sig_results.txt", sep="\t", sep2=c("", " ", ""))
save(gsea_df_sig, file= "/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/arm_week_ArmWeek/GSEA_wlzMAGs_sig_results.Rda")

more_tags_EdgeR_annotate <- left_join(more_tags_EdgeR, TableS5_mcseed_only_BP1)
more_tags_EdgeR_annotate <- left_join(more_tags_EdgeR_annotate, glossary)
write.csv(more_tags_EdgeR_annotate, "/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/arm_week_ArmWeek/all_tags_wlzMAGs.csv")
save(more_tags_EdgeR_annotate, file="/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/arm_week_ArmWeek/all_tags_wlzMAGs.Rda")

# Output contrast matrix:
write.csv(my.contrasts, "/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/arm_week_ArmWeek/contrast_matrix_modelArmWeekArmbyWeek.csv")
#################################################
####################################################

#Check which MAG and transcripts are driving carbohydrate enrichment 
lead_edg <- list()
seq <- seq(from=1, to=nrow(gsea_df_sig), by=1)
for(i in seq){
  name <- gsea_df_sig[i,c("pathway")]
  df <- cbind(gsea_df_sig[i,c(1:9,11:ncol(gsea_df_sig))],locus_tag=unlist(gsea_df_sig$leadingEdge[i]))
  df <- separate(df, locus_tag, c("locus","tag"), sep = "_", remove = FALSE)
  df <-  left_join(df, MAG_taxa[,c("MAG","locus","MAG_ID","Genus","Species") ])
  lead_edg[[name]] <- df # Get annotations for members of an enriched GSEA pathway
}

names(lead_edg) #Get a list of the enriched carbohydrate pathways and the annotations for their leading edges
rm(leading_edge_transcripts)
leading_edge_transcripts <- plyr::ldply(lead_edg, data.frame) #turn list into df with elements of the list stacked vertically. 
leading_edge_transcripts[,c(".id")] <- NULL #get rid of duplicate column
# Make a fig-friently MAG name
leading_edge_transcripts$MAG_ID_species <- paste(leading_edge_transcripts$MAG_ID,leading_edge_transcripts$Species)
save(leading_edge_transcripts, file="/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/arm_week_ArmWeek/gsea_leading_edges_all.Rda")
write.csv(leading_edge_transcripts, "/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/arm_week_ArmWeek/gsea_leading_edges_all.csv")

# Narrow the leading edge list down to just mcSEED 
leading_edge_transcripts_mcseed <- leading_edge_transcripts[leading_edge_transcripts$path=="mcSEED_phenotype",]
write.csv(leading_edge_transcripts_mcseed, "/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/arm_week_ArmWeek/gsea_leading_edges_mcSEED.csv")

# Filter the leading edges to carbs
leading_edge_transcripts_carbs <- leading_edge_transcripts_mcseed[leading_edge_transcripts_mcseed$Category=="carbohydrates" & !is.na(leading_edge_transcripts_mcseed$Category), ]
write.csv(leading_edge_transcripts_carbs, "/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/arm_week_ArmWeek/gsea_leading_edges_carbs.csv")

#For output
out <- leading_edge_transcripts
out[,c("leadingEdge")] <- NULL
colnames(out)[1] <- "group" #Rename the contrast so that this data can be merged with the annoations
#Merge with DGE data
out <- left_join(out, more_tags_EdgeR_annotate)
write.csv(out, "/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/arm_week_ArmWeek/gsea_leading_edges_all.csv")

#Annotate MDCF2 for pos NES and RUSF for neg NES
out$group_effect <- "hold"
out[out$NES>0,]$group_effect <- "MDCF2" #If NES > 1 then set to Q4
#There are no carb enrichements with neg NES; therefore the following line is not needed
out[out$NES<0,]$group_effect <- "RUSF" #If NES < 1 then set to Q1

# Count number of mcSEED pathways:
count_phenotype <- distinct(out, group, pathway, .keep_all = TRUE)
# Count Significant total:
dplyr::count(count_phenotype, Phenotype)
# Count Significant pathways for Q4/Q1:
dplyr::count(count_phenotype, Phenotype, group_effect)

# Limit the df to only the carbs
out


###################
#### Lollipop chart of MDCF2 vs RUSF carbohydrate pathway enrichement:
path_sig_select_carbs <- gsea_df_sig[gsea_df_sig$Category=="carbohydrates" & !is.na(gsea_df_sig$Category), ]
View(path_sig_select_carbs)

# Make -log10(P)
path_sig_select_carbs$log10pval <- -log10(path_sig_select_carbs$pval)
# Make sign -log10(P)
path_sig_select_carbs$Snlog10pval <- path_sig_select_carbs$log10pval
path_sig_select_carbs[path_sig_select_carbs$NES<0, c("Snlog10pval")] <- -(path_sig_select_carbs[path_sig_select_carbs$NES<0 ,c("log10pval")])
path_sig_select_carbs <-  path_sig_select_carbs[order(path_sig_select_carbs$padj, decreasing=FALSE),] #Order by padj
path_sig_select_carbs$Description <- factor(path_sig_select_carbs$Description, levels = path_sig_select_carbs$Description)


# B. Make chart
coefficient <- factor(path_sig_select_carbs$contrast)
pdf("/scratch/jglab/dmwebber/RNAseq/R_reanalysis_210826_onward/out/arm_week_ArmWeek/mcSEEDlevel_phenotype_CarbPathwaysEnrich.pdf", width = 5, height = 1.5)
ggplot(path_sig_select_carbs, aes(x = log10pval, y = Description, color = coefficient)) +
  geom_segment(aes(x = 0, y = Description, xend = log10pval, yend = Description), color = "grey") +
  geom_point(aes()) +
  scale_color_brewer(palette="Dark2") +
  theme_classic() +
  theme(axis.title.x = element_text(face="bold", size=7),
        axis.text.x  = element_text(size=6),
        axis.title.y = element_text(face="bold", size=7),
        axis.text.y  = element_text(size=6) +
          theme(legend.position = "right"))
dev.off()