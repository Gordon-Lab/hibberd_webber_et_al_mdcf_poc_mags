#rm(list=ls()) #Clear out old variables
# vignette: https://cran.r-project.org/web/packages/multilevelTools/vignettes/lmer-vignette.html
library(openxlsx) 
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
library(lme4)
library(lmerTest)
library(extraoperators)
library(JWileymisc)
library(multilevelTools)
library(caret)
library(R.utils)
library(Hmisc) #rcorr function
library(hrbrthemes)
library(viridis)
library(ggplot2)
library(ggpubr) #compare_means


# Bring in and summarize data
setwd("C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/code for deposition/hibberd_webber_et_al_mdcf_poc_mags/diet_intake")

datestring <- format(Sys.time(), "%y%m%d_%H%M")
study <- "MDCF_POC"
experiment <- "diet_intake"

cbbPalette <- c("#999999", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "darkred", "darkblue", "darkgrey", "white")

map.quartiles <- read.xlsx("/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/210810_1323_MDCF_POC_participant_beta_WLZ_quartiles.xlsx", sheet=1, startRow = 2, colNames = TRUE) %>% select(PID:beta_quartile)
map.fecal.timeline <- read.xlsx("/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/210223_MDCF_POC_fecal_plasma_metadata.xlsx", sheet=1, startRow = 1, colNames = TRUE)
metadata <- left_join(map.fecal.timeline, map.quartiles) #combine quartile and fecal sample timeline info 
metadata$sample_number <- as.numeric(metadata$sample_number)
Q1 <- metadata[metadata$beta_quartile=="Q1","SID"]
Q4 <- metadata[metadata$beta_quartile=="Q4","SID"]

metadata %>% 
  group_by(study_arm, beta_quartile) %>%
  summarize(mean(coef))

#Load GeTMM normalized counts: rpk_norm_counts_5_28_2022.Rdata
(load("C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/May2022_Fig4_follow_up_analysis/rpk_norm_counts_5_28_2022.Rdata"))
class(rpk_norm) 
#[1] "matrix" "array"
dim(rpk_norm)
#[1] 1929056     350

#Load fractional expression counts (GeTMM gene transcript counts / sum of all GeTMM transcript counts for MAG)
(load("C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/rpk_norm_fractional_expression_of_7749_filtered_genes.Rdata"))
# [1] "rpk_norm_WLZ_frac"

#Load RNA annotations
(load("C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/annotations_062621.Rda"))
dim(annotations)

# import supplement intake dataset
supp_intake <- read.xlsx("/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/Clinical_data/Daily Supplement Intake dataset.xlsx", sheet=1, startRow = 1, colNames = TRUE) %>% select(!(Group:IDF_Date))
supp_intake <- left_join(supp_intake, unique(select(metadata, PID, study_arm, beta_quartile)))

######load vsd mag abundance: Filters applied
#load MAG abundance
MAG_vst <- read.table("C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/MAGs/210614_0934_MDCF_POC_MAG_1000_counts_vst_filt.txt", sep= "\t", header = TRUE, row.names = 1)

#Load MAG data
MAG_meta_data <- read.xlsx("C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/reanalysis_210826_onward/input/210617_MDCF_POC_MAG_tables_rFriendly.xlsx", sheet=3)

#Import RNAseq data
tableS10 <- read.xlsx("C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/Volcano plot/TableS10_rIntake.xlsx", sheet=1)
WLZ_locus_tags <- tableS10 %>% filter(contrast == "arm") %>% select("locus_tag","MAG_ID","genus","species") %>% distinct(locus_tag)
locus_link <- tableS10 %>% filter(contrast == "arm") %>% select("locus_tag","MAG_ID","genus","species") %>% separate(locus_tag, c("locus","tag")) %>% distinct(locus, .keep_all = TRUE) %>% select(-c("tag"))

# import FFQ data
FFQ <- read.csv("FFQ.csv")
  #Replace NA with 0
  # FFQ[is.na(FFQ_raw)] <- 0
  # FFQ <- select_if(FFQ, is.numeric)

# import monosaccharide data
monos <- read.xlsx("/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/glycans/2021_10_4_Monosaccharide and Linkage analysis of POC study and Diets.xlsx", sheet=1, startRow = 2, colNames = TRUE)
  #Get rid of the extra column
  monos <- select(monos,!"Sample.number") 
  #Get rid of rows that have multiple sample replicates
  monos <- monos %>% filter(!grepl('_', Sample.ID))
  
# import linkage data
  linkage <- read.xlsx("/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/glycans/2021_10_4_Monosaccharide and Linkage analysis of POC study and Diets.xlsx", sheet=4, startRow = 2, colNames = TRUE, check.names = TRUE) %>% 
      select(-c("Sample.number")) %>%
      filter(!grepl('_', Sample.ID)) %>%
      column_to_rownames("Sample.ID")

# Load PUL locus_tag assignments
  PUL_tags <- read.csv("C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/PUL_locus_tags.csv")
  #Join these with the locus_link file
  
# Import MAG_link
  MAG_link <- read.csv("C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/MAG_link.csv") %>% distinct() 
  
# Import full set of MAG-locus links
  (load("C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/Input/MAG_locus_link.Rda"))
  dim(MAG_locus_link)
  #Join the PUL locus information onto the locus_link file
  PUL_locus_tags <- PUL_tags %>% separate(locus_tag,c("locus","tag"), sep="_", remove=FALSE) %>% left_join(MAG_locus_link, by="locus") %>% left_join(MAG_link, by="MAG") %>%
    mutate(locus_PUL= paste(locus, PUL, sep="_")) 
  #write.csv(PUL_locus_tags, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/PUL_locus_tags.csv")
    
######load vsd counts: Filters applied A) filtered out genes that were not part of 837 MAGs in Matt's abundance filtered set, 
  # B) filtered out genes annotated as rRNA 
  # C) set transcripts to zero for MAGs with abundance <.5 TPM 
  # D) Filtered out genes with <10 raw counts total among all samples
  (load(file="C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/May2022_Fig4_follow_up_analysis/210621_MDCF_POC_transcript_counts_vst_filtered_SumCtLT10_wTranscriptZeroing.Rda"))
  class(vsd)
  #[1] "DESeqTransform"
  #attr(,"package")
  #[1] "DESeq2"
  vsd <- assay(vsd)
  class(vsd)
  #[1] "matrix" "array"
  
# Exclude rows where either left1 or left2 is missing(NA)
supp_intake <- supp_intake[!(is.na(supp_intake$taken1) | is.na(supp_intake$taken2)),]
supp_intake$taken_daily <- supp_intake$taken1 + supp_intake$taken2

#Make a dummy variable
supp_intake$study_arm_num[supp_intake$study_arm=="MDCF2"] <- 2
supp_intake$study_arm_num[supp_intake$study_arm=="RUSF"] <- 1

#Examine the intra-class correlation coefficient or ICC
iccMixed(
  dv = "taken_daily",
  id = "PID",
  data = supp_intake)

#Check for outliers with person and between person 
tmp <- meanDecompose(taken_daily ~ PID, data = supp_intake)

  #Plot between person supplement consumed
  plot(testDistribution(tmp[["taken_daily by PID"]]$X,
                      extremevalues = "theoretical", ev.perc = .001),
     varlab = "Between person supplement consumed")
  
  #Plot within person supplement consumed
  plot(testDistribution(tmp[["taken_daily by residual"]]$X,
                        extremevalues = "theoretical", ev.perc = .001),
       varlab = "Within person supplement consumed")

  # test the model: taken_daily ~ day + study_arm_num + study_arm_num*day + (1|PID)
##############################################################################
  m <- lmerTest::lmer(taken_daily ~ Day + study_arm + study_arm:Day + (1|PID), data = supp_intake)
  
  #Run and plot model diagnostics
  md <- modelDiagnostics(m, ev.perc = .001)
  plot(md, ask = FALSE, ncol = 2, nrow = 3)
  
  #Evaluate outliers by subsetting to the random effect extreme values and viewing the top most values (head)
  rm(mvextreme)
  mvextreme <- subset(md[["extremeValues"]],
                      EffectType == "Random Effect PID : (Intercept)")
  head(mvextreme,50)
  unique(mvextreme$PID)
  # P09C666 (RUSF group)
  
  modelPerformance(m)
  # Note related to marginal vs conditional Cohen's F2 effect size: "A conditional treatment effect (conditional F2) is the average effect of treatment on the individual. 
  # A marginal treatment effect is the average effect of treatment on the population."
  # R2 measures of the variance accounted for by the fixed effects (marginal R2) and from the fixed and random effects combined (conditional R2)
  
  #see the results of individual variables:
  summary(m)

  # test the model: taken_daily ~ day + beta_quartile + beta_quartile*day + (1|PID)
############################################################################## 
  
  MDCF2supp_intake <- supp_intake[supp_intake$beta_quartile == c("Q1", "Q4") & supp_intake$study_arm == "MDCF2",]
  dim(MDCF2supp_intake)
  m2 <- lmerTest::lmer(taken_daily ~ Day + beta_quartile + beta_quartile:Day + (1|PID), data = MDCF2supp_intake)
  
  #Run and plot model diagnostics
  md <- modelDiagnostics(m2, ev.perc = .001)
  plot(md, ask = FALSE, ncol = 2, nrow = 3)
  
  #Evaluate outlines by subsetting to the random effect extreme values and viewing the top most values (head)
  rm(mvextreme)
  mvextreme <- subset(md[["extremeValues"]],
                      EffectType == "Random Effect PID : (Intercept)")
  head(mvextreme,50)
  unique(mvextreme$PID)

  
  modelPerformance(m2)
  # Note related to marginal vs conditional Cohen's F2 effect size: "A conditional treatment effect (conditional F2) is the average effect of treatment on the individual. 
  # A marginal treatment effect is the average effect of treatment on the population."
  # R2 measures of the variance accounted for by the fixed effects (marginal R2) and from the fixed and random effects combined (conditional R2)
  
  #see the results of individual variables:
  summary(m2)
  
  ####################################################
  ####################################################
  ## FFQ: Use the endpoint questions from FFQ as a classifier of upper- vs lower-quartile membership
  ##############
# Add metadata, especially the WLZ response quartile to FFQ
  FFQ <-  left_join(FFQ, metadata)
    
# Limit FFQ to a) timepoint 7, b) MDCF2, c) upper and lower quartiles
  FFQ_t7_MDCF2 <- FFQ[FFQ$sample_number==7 & FFQ$study_arm == "MDCF2",]
  FFQ_t7_MDCF2 <- FFQ_t7_MDCF2[!is.na(FFQ_t7_MDCF2$PID),]
  rownames(FFQ_t7_MDCF2) <- FFQ_t7_MDCF2$SID
  FFQ_t7_MDCF2Q1_4 <- FFQ_t7_MDCF2[FFQ_t7_MDCF2$beta_quartile %in% c("Q1", "Q4"),]
  FFQ_t7_MDCF2Q1_4$beta_quartile <- as.factor(FFQ_t7_MDCF2Q1_4$beta_quartile)
# Test for different intake by group
  compare_means(FFQ118A ~ beta_quartile, data = FFQ_t7_MDCF2, method = "t.test", paired = FALSE)
  
  # Make a composite measure of grains, roots and tubers (like fg1):
  FFQ_t7_MDCF2Q1_4_grains_roots <-  select(FFQ_t7_MDCF2Q1_4, FFQ114a, FFQ115a, FFQ116a)
  FFQ_t7_MDCF2Q1_4_grains_roots$grains_roots_srv <- rowSums(FFQ_t7_MDCF2Q1_4_grains_roots)
  FFQ_t7_MDCF2Q1_4_grains_roots_bi <- select(FFQ_t7_MDCF2Q1_4, FFQ114, FFQ115,FFQ116)
  FFQ_t7_MDCF2Q1_4_grains_roots_bi$grains_roots_bi <- rowSums(FFQ_t7_MDCF2Q1_4_grains_roots_bi)
  # Make a composite measure of dairy products (fg3):
  FFQ_t7_MDCF2Q1_4_dairy <-  select(FFQ_t7_MDCF2Q1_4, FFQ102, FFQ103, FFQ105, FFQ106,FFQ125a)
  FFQ_t7_MDCF2Q1_4_dairy$dairy_srv <- rowSums(FFQ_t7_MDCF2Q1_4_dairy)
  FFQ_t7_MDCF2Q1_4_dairy_bi <-  select(FFQ_t7_MDCF2Q1_4, FFQ101, FFQ104, FFQ125, FFQ106,FFQ125a)
  FFQ_t7_MDCF2Q1_4_dairy_bi$dairy_bi <- rowSums(FFQ_t7_MDCF2Q1_4_dairy_bi)
  # Make a composite measure of fruits and veggies:
  FFQ_t7_MDCF2Q1_4_fruits_veg <-  select(FFQ_t7_MDCF2Q1_4, FFQ117a, FFQ119a, FFQ120a)
  FFQ_t7_MDCF2Q1_4_fruits_veg$fruits_veg_srv <- rowSums(FFQ_t7_MDCF2Q1_4_fruits_veg)
  FFQ_t7_MDCF2Q1_4_fruits_veg_bi <-  select(FFQ_t7_MDCF2Q1_4, FFQ117, FFQ119, FFQ120)
  FFQ_t7_MDCF2Q1_4_fruits_veg_bi$fruits_veg_bi <- rowSums(FFQ_t7_MDCF2Q1_4_fruits_veg_bi)
  
  # Separate the df into metadata and numeric data  
  FFQ_t7_MDCF2Q1_4_metadata <- FFQ_t7_MDCF2Q1_4 %>% select(sample_number:day, enrollment_age:beta_quartile)
  FFQ_t7_MDCF2Q1_4_FFQ <-  select(FFQ_t7_MDCF2Q1_4, !colnames(FFQ_t7_MDCF2Q1_4_metadata))
    #Check the SD to see if there is some level of variation in each predictor var
    FFQ_SD <-apply(FFQ_t7_MDCF2Q1_4_FFQ, 2, sd)
    # FFQ_mean <- apply(FFQ_t7_MDCF2Q1_4_FFQ, 2, mean)
    # FFQ_CV <- FFQ_SD/FFQ_mean
    FFQ_nonZero_SD <- FFQ_SD[FFQ_SD>0]
    #Limit the df to the variables with SD>0
    FFQ_t7_MDCF2Q1_4_FFQ <- FFQ_t7_MDCF2Q1_4_FFQ[,colnames(FFQ_t7_MDCF2Q1_4_FFQ) %in% names(FFQ_nonZero_SD)]
  FFQ_t7_MDCF2Q1_4_FFQ_qrt <- cbind(FFQ_t7_MDCF2Q1_4_FFQ,beta_quartile = FFQ_t7_MDCF2Q1_4_metadata$beta_quartile, 
                                    grains_roots_srv =FFQ_t7_MDCF2Q1_4_grains_roots$grains_roots_srv, 
                                    dairy_srv = FFQ_t7_MDCF2Q1_4_dairy$dairy_srv, fruits_veg_srv = FFQ_t7_MDCF2Q1_4_fruits_veg$fruits_veg_srv,
                                    grains_roots_bi = FFQ_t7_MDCF2Q1_4_grains_roots_bi$grains_roots_bi, dairy_bi = FFQ_t7_MDCF2Q1_4_dairy_bi$dairy_bi,
                                    fruits_veg_bi = FFQ_t7_MDCF2Q1_4_fruits_veg_bi$fruits_veg_bi)
  FFQ_t7_MDCF2Q1_4_FFQ_qrt$Sample.ID <- rownames(FFQ_t7_MDCF2Q1_4_FFQ_qrt)
  
  ######## ML modeling: MDCF2-WLZ quartile ~ FFQ_t7
  #### 1. Set random seed
  set.seed(100)
  
  #### 2.Create data partition
  trainRowNumbers <- createDataPartition(FFQ_t7_MDCF2Q1_4$beta_quartile, p=0.7, list=FALSE)
  
  #### 3. Partition the data
  # Create the training  dataset
  train_data <- FFQ_t7_MDCF2Q1_4_FFQ[trainRowNumbers,]
  train_metadata <- FFQ_t7_MDCF2Q1_4_metadata[trainRowNumbers,]
  train_data_y <- FFQ_t7_MDCF2Q1_4_FFQ_qrt[trainRowNumbers,]

    #Check the "training" fraction
    nrow(train_data)/nrow(FFQ_t7_MDCF2Q1_4_metadata)
    nrow(train_metadata)/nrow(FFQ_t7_MDCF2Q1_4_metadata)  
  
  # Create the test dataset
  test_data <- FFQ_t7_MDCF2Q1_4_FFQ[-trainRowNumbers,]
  test_metadata <- FFQ_t7_MDCF2Q1_4_metadata[-trainRowNumbers,]
  test_data_y <- FFQ_t7_MDCF2Q1_4_FFQ_qrt[-trainRowNumbers,]
  
    #Check the "test" fraction
    nrow(test_data)/nrow(FFQ_t7_MDCF2Q1_4_metadata)
    nrow(test_metadata)/nrow(FFQ_t7_MDCF2Q1_4_metadata)
    
    ############################ b. Apply a elastic net model
    # 1. Choose how the resampling is done to identify best parameters
    ctrl <- trainControl(
      # Choose the method
      method = 'repeatedcv',
      # Choose the number of folds
      number = 5,  #number of K-fold cross-validation. Try 10
      # Choose the number of times to repeat the cv
      repeats = 5,  #number of times to repeat K-fold cross validation. Try 10
      classProbs = TRUE)
      #Subsampling is only used for binomial sampling.    
    
      # 2. Train the ridge regression model, e.g. alpha = 0. for glmnet
      glmnetFit <- caret::train(
      # Set the formula
      beta_quartile ~ .,
      # Specify the data
      data = train_data_y, 
      # Choose the resampling method
      trControl = ctrl,
      # Choose the method (ada, bag, bagEarth, bagFDA, blackboost, cforest, ctree, ctree2, rf)
      method = 'gbm',
      # Choose the preprocessing, although this is extra
      family = "binomial", #other options: "gaussian", "binomial", "poisson", "multinomial", "cox", "mgaussian"
      # type.measure = "auc", #options: deviance (uses actual deviance), mae (uses mean absolute error),class (gives misclassification error), auc (for two-class logistic regression ONLY, gives area under the ROC curve)
      preProc = c('center', 'scale'),
      # Choose the number of parameters. For PLS, this is the number of PCs to choose
      tuneLength = 20, # try 50
      # tuneGrid = expand.grid(alpha = seq(0,1,length = 10), lambda = seq(0, 1, length = 10)), # tuneGrid = expand.grid(alpha = seq(0,.2,.01), lambda = seq(0.0001, 1, length = 1000))
      # metric = "Kappa" #Kappa, Accuracy
      )
      
    #Show statistics for the top 10 performing ML settings
    head(as.data.frame(glmnetFit$results) %>% arrange(desc(Kappa)),10)
    max(glmnetFit[["results"]][["Accuracy"]])
    
    # 3. Evaluate the crossvalidation
    ### i. Take a look at the model
    plot(glmnetFit)
    write.csv(as.data.frame(glmnetFit$results), "./R_output/anthro_disc_ML2_performanceA.csv")
    
    ### ii. Test the model
    glmnetPredict <- predict(glmnetFit, newdata = test_data_y)
    postResample(pred = glmnetPredict, obs=(as.factor(test_data_y$beta_quartile)))
    ML_performance <- confusionMatrix(data = glmnetPredict, reference = (as.factor(test_data_y$beta_quartile)))
    ML_performance
    write.csv(ML_performance$byClass, "ML_performanceB.csv")
    
    importance <- varImp(glmnetFit, scale=FALSE) #varImp is scaled to have a maximum value of 100, unless the scale argument of varImp.train is set to FALSE.
    importance <- importance$importance
    importance <- rownames_to_column(importance, var="SeqId") 

##############################################################################################
#####################
  # Check for a correlation between FFQ questions and fecal glycans in upper and lower quartile
  FFQ_monos_MDCF2_Q1_4 <- left_join(FFQ_t7_MDCF2Q1_4_FFQ_qrt, monos)
  rownames(FFQ_monos_MDCF2_Q1_4) <- FFQ_monos_MDCF2_Q1_4$Sample.ID
  FFQ_monos_MDCF2_Q1_4$Sample.ID <- NULL
  FFQ_monos_MDCF2_Q1_4$beta_quartile <- NULL
  corr_FFQ_monos_MDCF2_Q1_4 <- rcorr(as.matrix(FFQ_monos_MDCF2_Q1_4), type = "spearman") 
  r_FFQ_monos_MDCF2_Q1_4 <- corr_FFQ_monos_MDCF2_Q1_4[["r"]]
  p_FFQ_monos_MDCF2_Q1_4 <- corr_FFQ_monos_MDCF2_Q1_4[["P"]]
  r_FFQ_monos_MDCF2_Q1_4 <- r_FFQ_monos_MDCF2_Q1_4[rownames(r_FFQ_monos_MDCF2_Q1_4) %in% colnames(monos),colnames(r_FFQ_monos_MDCF2_Q1_4) %in% colnames(FFQ_t7_MDCF2Q1_4_FFQ_qrt) ]
  p_FFQ_monos_MDCF2_Q1_4 <- p_FFQ_monos_MDCF2_Q1_4[rownames(p_FFQ_monos_MDCF2_Q1_4) %in% colnames(monos),colnames(p_FFQ_monos_MDCF2_Q1_4) %in% colnames(FFQ_t7_MDCF2Q1_4_FFQ_qrt) ]
  
  # Filter the df to fg2, fg3, fg4
  p_FFQ_monos_MDCF2_Q1_4_fg <- p_FFQ_monos_MDCF2_Q1_4[,c("fg2","fg3","fg4","fg5", "fg6", "fg7")]
  r_FFQ_monos_MDCF2_Q1_4_fg <- r_FFQ_monos_MDCF2_Q1_4[,c("fg2","fg3","fg4","fg5", "fg6", "fg7")]
  
  #Apply p.adjust to fg2:fg4
  p_FFQ_monos_MDCF2_Q1_4_fg %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(p_FFQ_monos_MDCF2_Q1_4_fg))

  write.csv(r_FFQ_monos_MDCF2_Q1_4, "r_corr_FFQ_monos_MDCF2_Q1_4_spearman.csv")
  write.csv(p_FFQ_monos_MDCF2_Q1_4, "p_corr_FFQ_monos_MDCF2_Q1_4_spearman.csv")
  
  ## Make a dotplot correlating various monosaccharides with fg2 intake
  # Basic scatter plot
  ggplot(FFQ_monos_MDCF2_Q1_4, aes(x=wt, y=mpg)) + geom_point()
  
  FFQ_monos_MDCF2_Q1_4$fg2 <- as.factor(FFQ_monos_MDCF2_Q1_4$fg2)
  p1<-ggplot(FFQ_monos_MDCF2_Q1_4, aes(x=FFQ101, y=Galactose)) + geom_point() 
  p1
  
    ## Spot check results from correlation analysis by using T-test
    ## Perform T-test of mon ~ FFQ
    #fg2
    #First convert from wide to tall df
    FFQ_monos_MDCF2_Q1_4$fg2 <- as.factor(FFQ_monos_MDCF2_Q1_4$fg2)
    FFQ_monos_MDCF2_Q1_4_tall <- select(FFQ_monos_MDCF2_Q1_4, "fg2", Glucose:Ribose) %>% pivot_longer(!fg2, names_to="glycan", values_to = "quant")
    #Run T-test
    compare_means(quant ~ fg2, data = FFQ_monos_MDCF2_Q1_4_tall, group.by = "glycan", method = "t.test", paired = FALSE)

  #####################
  # Check for a correlation between FFQ questions and fecal linkages in upper and lower quartile
    FFQ_linkage_MDCF2_Q1_4 <- left_join(FFQ_t7_MDCF2Q1_4_FFQ_qrt, linkage)
    rownames(FFQ_linkage_MDCF2_Q1_4) <- FFQ_linkage_MDCF2_Q1_4$Sample.ID
    FFQ_linkage_MDCF2_Q1_4$Sample.ID <- NULL
    FFQ_linkage_MDCF2_Q1_4$beta_quartile <- NULL
    corr_FFQ_linkage_MDCF2_Q1_4 <- rcorr(as.matrix(FFQ_linkage_MDCF2_Q1_4), type = "spearman") 
    r_FFQ_linkage_MDCF2_Q1_4 <- corr_FFQ_linkage_MDCF2_Q1_4[["r"]]
    p_FFQ_linkage_MDCF2_Q1_4 <- corr_FFQ_linkage_MDCF2_Q1_4[["P"]]
    r_FFQ_linkage_MDCF2_Q1_4 <- r_FFQ_linkage_MDCF2_Q1_4[rownames(r_FFQ_linkage_MDCF2_Q1_4) %in% colnames(linkage),colnames(r_FFQ_linkage_MDCF2_Q1_4) %in% colnames(FFQ_t7_MDCF2Q1_4_FFQ_qrt) ]
    p_FFQ_linkage_MDCF2_Q1_4 <- p_FFQ_linkage_MDCF2_Q1_4[rownames(p_FFQ_linkage_MDCF2_Q1_4) %in% colnames(linkage),colnames(p_FFQ_linkage_MDCF2_Q1_4) %in% colnames(FFQ_t7_MDCF2Q1_4_FFQ_qrt) ]
    
  #Limit linkages to those 14 that are in Fig4
    Sig14linkages <- c("T-F-Arabinose", "2-F-Arabinose", "5-F-Arabinose", "2,3-F-Arabinose", "T-P-Xylose", "2-Xylose", "T-GlcA", "T-GalA", "6-Galactose", "2,4,6-Glucose", "3-Mannose", "4-Mannose", "4,6-Mannose", "3,4,6-Mannose")
    r_FFQ_14linkage_MDCF2_Q1_4 <- r_FFQ_linkage_MDCF2_Q1_4[Sig14linkages, colnames(r_FFQ_linkage_MDCF2_Q1_4) %in% colnames(FFQ_t7_MDCF2Q1_4_FFQ_qrt) ]
    p_FFQ_14linkage_MDCF2_Q1_4 <- p_FFQ_linkage_MDCF2_Q1_4[Sig14linkages, colnames(p_FFQ_linkage_MDCF2_Q1_4) %in% colnames(FFQ_t7_MDCF2Q1_4_FFQ_qrt) ]
    
    write.csv(r_FFQ_14linkage_MDCF2_Q1_4, "r_corr_FFQ_linkages_MDCF2_Q1_4_spearman.csv")
    write.csv(p_FFQ_14linkage_MDCF2_Q1_4, "p_corr_FFQ_linkages_MDCF2_Q1_4_spearman.csv")
    
    # Filter the df to fg2, fg3, fg4
    p_FFQ_linkages_MDCF2_Q1_4_fg <- p_FFQ_14linkage_MDCF2_Q1_4[,c("fg2","fg3","fg4","fg5", "fg6", "fg7")]
    r_FFQ_linkages_MDCF2_Q1_4_fg <- r_FFQ_14linkage_MDCF2_Q1_4[,c("fg2","fg3","fg4","fg5", "fg6", "fg7")]
    
    #Apply p.adjust to fg2:fg4
    pAdj_FFQ_linkages_MDCF2_Q1_4_fg <- p_FFQ_linkages_MDCF2_Q1_4_fg %>% as.vector %>% p.adjust(method='fdr') %>% matrix(ncol=ncol(p_FFQ_linkages_MDCF2_Q1_4_fg))
    rownames(pAdj_FFQ_linkages_MDCF2_Q1_4_fg) <- rownames(p_FFQ_linkages_MDCF2_Q1_4_fg)
    colnames(pAdj_FFQ_linkages_MDCF2_Q1_4_fg) <- colnames(p_FFQ_linkages_MDCF2_Q1_4_fg)
    
    write.csv(r_FFQ_linkages_MDCF2_Q1_4_fg, "r_corr_FFQ_linkages_MDCF2_Q1_4_spearman_fg.csv")
    write.csv(pAdj_FFQ_linkages_MDCF2_Q1_4_fg, "pAdj_corr_FFQ_linkages_MDCF2_Q1_4_spearman_fg.csv")
    
    
  # Check for a correlation between FFQ questions and fecal glycans in all participants at T7
  ############################################  
    # First, Get rid of non-numeric columns
    FFQ_t7_MDCF2 <- select_if(FFQ_t7_MDCF2,is.numeric)
    #Check the SD to see if there is some level of variation in each predictor var
    FFQ_SD <-sapply(FFQ_t7_MDCF2, sd)
    FFQ_SD <- FFQ_SD[!is.na(FFQ_SD)]
    FFQ_nonZero_SD <- FFQ_SD[FFQ_SD>0]
    #Limit the df to the variables with SD>0
    FFQ_t7_MDCF2 <- FFQ_t7_MDCF2[,colnames(FFQ_t7_MDCF2) %in% names(FFQ_nonZero_SD)]
    FFQ_t7_MDCF2$Sample.ID <- rownames(FFQ_t7_MDCF2)
    
    # Include only columns from the, "how many times questions"
    # FFQ_t7_MDCF2 <- data.frame(FFQ_t7_MDCF2[,grepl("a", colnames(FFQ_t7_MDCF2),)],Sample.ID = FFQ_t7_MDCF2$Sample.ID)
    
    #Run the actual correlation for fg2-fg7
    FFQ_monos_MDCF2 <- left_join(FFQ_t7_MDCF2, monos)
    rownames(FFQ_monos_MDCF2) <- FFQ_monos_MDCF2$Sample.ID
    FFQ_monos_MDCF2$Sample.ID <- NULL
    corr_FFQ_monos_MDCF2 <- rcorr(as.matrix(FFQ_monos_MDCF2), type = "spearman") 
    r_FFQ_monos_MDCF2 <- corr_FFQ_monos_MDCF2[["r"]]
    p_FFQ_monos_MDCF2 <- corr_FFQ_monos_MDCF2[["P"]]

    r_FFQ_monos_MDCF2 <- r_FFQ_monos_MDCF2[rownames(r_FFQ_monos_MDCF2) %in% colnames(monos),colnames(r_FFQ_monos_MDCF2) %in% colnames(FFQ_t7_MDCF2Q1_4_FFQ_qrt) ]
    p_FFQ_monos_MDCF2 <- p_FFQ_monos_MDCF2[rownames(p_FFQ_monos_MDCF2) %in% colnames(monos),colnames(p_FFQ_monos_MDCF2) %in% colnames(FFQ_t7_MDCF2Q1_4_FFQ_qrt) ]
    write.csv(r_FFQ_monos_MDCF2, "r_corr_FFQ_monos_MDCF2.csv")
    write.csv(p_FFQ_monos_MDCF2, "p_corr_FFQ_monos_MDCF2.csv")
  
     ###############
   ## Make a dotplot showing Food group consumption
    FFQ_t7_MDCF2Q1_4_FFQ_qrt$fg2 <- as.factor(FFQ_t7_MDCF2Q1_4_FFQ_qrt$fg2)
    p <- ggplot(data=FFQ_t7_MDCF2Q1_4_FFQ_qrt, aes(x=fg2, group=beta_quartile, fill=beta_quartile)) +
      geom_density(adjust=.5, alpha=.4) +
      theme_ipsum()
    p
  
  ## Make a dotplot correlating various monosaccharides with fg2 intake
    FFQ_monos_MDCF2$fg2 <- as.factor(FFQ_monos_MDCF2$fg2)
        p2<-ggplot(FFQ_monos_MDCF2, aes(x=fg2, y=GalA)) + 
      geom_dotplot(binaxis='y', stackdir='center')
    p2
 
    ## Perform T-test of mon ~ FFQ
    #First convert from wide to tall df
    FFQ_monos_MDCF2$fg2 <- as.factor(FFQ_monos_MDCF2$fg2)
    FFQ_monos_MDCF2_tall <- select(FFQ_monos_MDCF2, "fg2", Glucose:Ribose) %>% pivot_longer(!fg2, names_to="glycan", values_to = "quant")
    #Run T-test
    compare_means(quant ~ fg2, data = FFQ_monos_MDCF2_tall, group.by = "glycan", method = "t.test", paired = FALSE)
    
########################################################################################
############### Run linear mixed effects model to determine the effect of linkage availability on MAG abundance: 
    ########### MAG VST abundance ~ linkage abundance + time + (linkage x time) + PID
    ####Inputs:
    # linkage data
    dim(linkage)
    #[1] 179  50
    head(linkage,n=2)
    linkage <- column_to_rownames(linkage,"Sample.ID")
    #    Sample.ID T-Glucose 4-Glucose 6-Glucose 3-Glucose/3-Galactose 2-Glucose 4,6-Glucose 3,4-Glucose 2,4-Glucose 3,4,6-Glucose 2,4,6-Glucose T-Galactose
    #1 P01C6661102   1062908   8009238  42694.51  406194.1              37934.24   187777.37    47588.80    60024.18      13445.56      1104.102   2078927.0
    
    dim(MAG_vst)
    # [1] 837 707
    head(MAG_vst)
    #                        P01C5551102 P01C5551107 P01C5551108 P01C5551104 P01C5551105 P01C5551106 P01C6661102 P01C6661107 P01C6661108 P01C6661104 P01C6661105
    #P01C555_maxbin.024.mod2    12.65557    11.59686    8.966784    17.16807    18.41875    15.04337    7.649803     6.66174    8.769358    9.789561    10.31465    
    
    head(metadata)
    #          SID     PID gender enrollment_age baseline_age study_arm age_at_sampling study_week sample_number study_phase  study_phase_alt sample_type
    #1 P01C5551101 P01C555 Female        17.3589     17.58904      RUSF        17.39178         -1             1      run_in           run_in       fecal
    
    head(MAG_meta_data)
    #                       MAG   MAG_ID        Genus                 Species         beta      lower_CI     upper_CI         SEM      p_value
    #1  P01C555_maxbin.024.mod2 MAG_0619 Anaerostipes     Anaerostipes_hadrus  0.009022589 -4.471471e-03  0.022496841 0.006880725 0.1907348258
    
    # Filter MAG_vst to WLZ-associated MAGs
    MAG_vst_WLZassoc <- MAG_vst %>% rownames_to_column("MAG") %>%
      left_join(MAG_meta_data) %>%
      filter(MAG_abd_FDR_adjusted < 0.05) %>%
      select(-c("beta","lower_CI","upper_CI","SEM","p_value","MAG_abd_FDR_adjusted"))
    
    # transpose MAG_vst_WLZassoc so that rownames=SID
    MAG_vst_linkage <- MAG_vst_WLZassoc %>%
      select(-c("MAG","Genus","Species")) %>%
      column_to_rownames("MAG_ID")
    
    MAG_vst_linkage_t <- MAG_vst_linkage %>% t() %>% as.data.frame %>% rownames_to_column("SID") %>%
      inner_join(linkage %>% rownames_to_column("SID")) %>%
      left_join(metadata %>% select("SID","PID","study_week")) %>% 
      mutate(across(-c("SID","PID","study_week"), scale)) #Scale the df
    
    #write.csv(MAG_vst_linkage_t, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/linkage_MAG_abdn_df_220818.csv")
    MAG_vst_linkage_t <- read.csv("C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/linkage_MAG_abundance_df_220818.csv") %>% select(-c("X"))
    
    #MAG_vst_linkage_t$study_week <- as.factor(MAG_vst_linkage_t$study_week)
    ##Plot
    ggplot(data  = MAG_vst_linkage_t,
           aes(x = MAG_0073,
               y = X.Hexose,
               color = study_week))+
      geom_point()+
      theme_minimal()+
      #geom_line(aes(group = PID)) +
      labs(title = "MAG_0073 abundances vs. Hexose")
    
    ###############################
    tmp <- meanDecompose(MAG_0073 ~ PID, data = MAG_vst_linkage_t)
    
    #Plot between person "MAG_0018 by X5.F.Arabinose
    plot(testDistribution(tmp[["MAG_0073 by PID"]]$X,
                          extremevalues = "theoretical", ev.perc = .001),
         varlab = "Between person MAG_0073 abundance")
    
    #Plot within person supplement consumed
    plot(testDistribution(tmp[["MAG_0073 by residual"]]$X,
                          extremevalues = "theoretical", ev.perc = .001),
         varlab = "Within person MAG_0073 abundance")
    
    # test the model: MAG VST ~ linkage quant + week + linkage quant:week + (1|PID)
    ##############################################################################
    #Where i represents MAGs and j represents linkages
    i_length=length(rownames(MAG_vst_linkage))
    j_length=length(colnames(linkage))
    i=1
    j=1
    lmer_bulk <- data.frame()
    for (i in seq(1,i_length,1)) 
    {
    MAG_name=rownames(MAG_vst_linkage)[i]
      for (j in seq(1,j_length,1)) 
      {
      linkage_name=colnames(linkage)[j]
      f1 <-paste(MAG_name, "~", linkage_name, "+", "study_week", "+", linkage_name, ":study_week", "+", "(1|PID)", sep=" ")
      m <- lmerTest::lmer(f1, data = MAG_vst_linkage_t)
      m_coef <- summary(m)[["coefficients"]] %>% as.data.frame() %>% rownames_to_column("coef") %>% select(-c("coef")) %>% mutate(MAG = MAG_name, linkage = linkage_name) %>% cbind(coef = c("intercept","beta1", "beta2", "beta3"))
      lmer_bulk <- lmer_bulk %>% bind_rows(m_coef)
      }
    }
    lmer_bulk <- lmer_bulk %>% rename("Pr(>|t|)"="p_value")
    lmer_bulk_1 <- lmer_bulk %>% left_join(select(MAG_meta_data, -c("lower_CI","upper_CI","SEM","p_value")), by = c("MAG" = "MAG_ID")) %>% rename("beta"="MAG_abd_beta") 
    lmer_bulk_linkageXweek <- lmer_bulk_1 %>% filter(coef=="beta3") %>% mutate(q_value = p.adjust(p_value))
    lmer_bulk_linkage <- lmer_bulk_1 %>% filter(coef=="beta1") %>% mutate(q_value = p.adjust(p_value))
    write.csv(lmer_bulk_linkageXweek, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/mag_response_to_linkage_linkageXweek_week_beta3_220818.csv")
    write.csv(lmer_bulk_linkage, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/mag_response_to_linkage_linkageXweek_week_beta1_220818.csv")
    
    # test the model: MAG VST ~ linkage quant + (1|PID)
    ##############################################################################
    #Where i represents MAGs and j represents linkages
    rm(i_length,j_length,f2,m2,m_coef2)
    i_length=length(rownames(MAG_vst_linkage))
    j_length=length(colnames(linkage))
    i=1
    j=1
    lmer_bulk2 <- data.frame()
    for (i in seq(1,i_length,1)) 
    {
      MAG_name=rownames(MAG_vst_linkage)[i]
      for (j in seq(1,j_length,1)) 
      {
        linkage_name=colnames(linkage)[j]
        f2 <-paste(MAG_name, "~", linkage_name, "+", "(1|PID)", sep=" ")
        m2 <- lmerTest::lmer(f2, data = MAG_vst_linkage_t)
        m_coef2 <- summary(m2)[["coefficients"]] %>% as.data.frame() %>% rownames_to_column("coef") %>% select(-c("coef")) %>% mutate(MAG = MAG_name, linkage = linkage_name) %>% cbind(coef = c("intercept","beta1"))
        lmer_bulk2 <- lmer_bulk2 %>% bind_rows(m_coef2)
      }
    }
    lmer_bulk2 <- lmer_bulk2 %>% rename("Pr(>|t|)"="p_value")
    lmer_bulk2_linkage <- lmer_bulk2 %>% filter(coef=="beta1") %>% mutate(q_value = p.adjust(p_value))
    dim(lmer_bulk2_linkage)
    write.csv(lmer_bulk2_linkage, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/mag_response_to_linkage_220818.csv")
    
    # test the model: MAG VST ~ linkage quant + week + (1|PID)
    ##############################################################################
    #Where i represents MAGs and j represents linkages
    i_length=length(rownames(MAG_vst_linkage))
    j_length=length(colnames(linkage))
    i=1
    j=1
    lmer_bulk3 <- data.frame()
    for (i in seq(1,i_length,1)) 
    {
      MAG_name=rownames(MAG_vst_linkage)[i]
      for (j in seq(1,j_length,1)) 
      {
        linkage_name=colnames(linkage)[j]
        f3 <-paste(MAG_name, "~", linkage_name, "+", "study_week", "+", "(1|PID)", sep=" ")
        m3 <- lmerTest::lmer(f3, data = MAG_vst_linkage_t)
        m_coef3 <- summary(m3)[["coefficients"]] %>% as.data.frame() %>% rownames_to_column("coef") %>% select(-c("coef")) %>% mutate(MAG = MAG_name, linkage = linkage_name) %>% cbind(coef = c("intercept","beta1", "beta2"))
        lmer_bulk3 <- lmer_bulk3 %>% bind_rows(m_coef3)
      }
    }
    lmer_bulk3 <- lmer_bulk3 %>% rename("Pr(>|t|)"="p_value")
    lmer_bulk3_linkage <- lmer_bulk3 %>% filter(coef=="beta1") %>% mutate(q_value = p.adjust(p_value))
    lmer_bulk3_linkage <- lmer_bulk3_linkage %>% left_join(select(MAG_meta_data, -c("lower_CI","upper_CI","SEM","p_value")), by = c("MAG" = "MAG_ID")) %>% rename("beta"="MAG_abd_beta") 
    dim(lmer_bulk3_linkage)
    write.csv(lmer_bulk3_linkage, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/mag_response_to_linkage_week_beta1_220818.csv")
 
##############
## model the effects of linkages on transcripts
# A. Filter vst counts to those that pass DGE filters and are from WLZ-associated MAGs.
    # Make a list of the significant linkages from Table 15a
    T15a_linkage_names <- c("X4.6.Mannose", "T.GalA", "T.GlcA", "T.P.Xylose", "X2.4.6.Glucose", "T.F.Arabinose", "X4.Mannose", "X5.F.Arabinose", "X6.Galactose", "X2.3.F.Arabinose", "X2.Xylose", "X3.4.6.Mannose", "X3.Mannose", "X2.F.Arabinose")
    T15a_linkage_NOnames <- colnames(linkage[,!colnames(linkage) %in% T15a_linkage_names])
    
    vsd_wlz <-  vsd[rownames(vsd) %in% WLZ_locus_tags$locus_tag,]
    # transpose df
    MAG_transcript_linkage_t <- vsd_wlz %>% t() %>% as.data.frame %>% rownames_to_column("SID") %>% inner_join(MAG_vst_linkage_t) %>% select(-c(T15a_linkage_NOnames))

    # test the model: transcript VST ~ linkage quant + week + linkage quant:week + (1|PID)
    ##############################################################################
    #Where i represents genes/transcripts and j represents linkages
    i_length=length(rownames(vsd_wlz))
    j_length=length(T15a_linkage_names)
    i=1
    j=1
    lmer_bulk <- data.frame()
    for (i in seq(1,i_length,1)) 
    {
      transcript_name=rownames(vsd_wlz)[i]
      for (j in seq(1,j_length,1)) 
      {
        linkage_name=T15a_linkage_names[j]
        f1 <-paste(transcript_name, "~", linkage_name, "+", "study_week", "+", linkage_name, ":study_week", "+", "(1|PID)", sep=" ")
        m <- lmerTest::lmer(f1, data = MAG_transcript_linkage_t)
        m_coef <- summary(m)[["coefficients"]] %>% as.data.frame() %>% rownames_to_column("coef") %>% select(-c("coef")) %>% mutate(transcript = transcript_name, linkage = linkage_name) %>% cbind(coef = c("intercept","beta1", "beta2", "beta3"))
        lmer_bulk <- lmer_bulk %>% bind_rows(m_coef)
      }
    }
    lmer_bulk <- lmer_bulk %>% rename("Pr(>|t|)"="p_value", "transcript"="tag")
    lmer_bulk <- lmer_bulk %>% rename("transcript"="locus_tag")
    write.csv(lmer_bulk, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/transcript_response_to_linkage_linkageXweek_week_220819.csv")
    
      lmer_bulk_1 <- lmer_bulk %>% left_join(annotations)
      lmer_bulk_linkageXweek <- lmer_bulk_1 %>% filter(coef=="beta3") %>% mutate(q_value = p.adjust(p_value))
      lmer_bulk_linkage <- lmer_bulk_1 %>% filter(coef=="beta1") %>% mutate(q_value = p.adjust(p_value))
      write.csv(lmer_bulk_linkageXweek, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/transcript_response_to_linkage_linkageXweek_week_beta3_220818.csv")
      write.csv(lmer_bulk_linkage, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/transcript_response_to_linkage_week_beta1_220818.csv")

#########################################################
######### model the effects of linkages on GeTMM factional expression
      head(rpk_norm_WLZ_frac)
      #                    P01C5551102 P01C5551105 P01C5551107 P01C6661102 P01C6661105 P01C6661107 P01C7771102 P01C7771105 P01C7771107 P01C8881102 P01C8881105 P01C8881107
      # MBIN0035D7B6_05620  0.06826254   0.5585589 0.232495356  0.33946184  0.03701932   0.7227685  0.07407188  0.05909928  0.02457301   0.8678454   0.5259785  0.01021300
      # MBIN0035D7B6_10265  0.01693486   0.0000000 0.008239778  0.01203074  0.02371296   0.0000000  0.01400082  0.02324911  0.02276735   0.0000000   0.0093205  0.01286333
      
      rpk_norm_WLZ_frac_wlz <-  rpk_norm_WLZ_frac[rownames(rpk_norm_WLZ_frac) %in% WLZ_locus_tags$locus_tag,]
      # transpose df, join with linkage data, get rid of all linkages except for those in Fig15a, center and scale the df
      rpk_norm_WLZ_frac_wlz_t <- rpk_norm_WLZ_frac_wlz %>% t() %>% as.data.frame %>% rownames_to_column("SID") %>% inner_join(MAG_vst_linkage_t) %>% select(-c(T15a_linkage_NOnames)) %>% mutate(across(-c("SID","PID","study_week"), scale))

      # test the model: transcript VST ~ linkage quant + week + linkage quant:week + (1|PID)
      ##############################################################################
      #Where i represents genes/transcripts and j represents linkages
      rm(i_length, j_length, transcript_name, linkage_name, f1, m, m_coef, lmer_bulk)
      i_length=length(rownames(rpk_norm_WLZ_frac_wlz))
      j_length=length(T15a_linkage_names)
      i=1
      j=1
      lmer_bulk <- data.frame()
      for (i in seq(1,i_length,1)) 
      {
        transcript_name=rownames(rpk_norm_WLZ_frac_wlz)[i]
        for (j in seq(1,j_length,1)) 
        {
          linkage_name=T15a_linkage_names[j]
          f1 <-paste(transcript_name, "~", linkage_name, "+", "study_week", "+", linkage_name, ":study_week", "+", "(1|PID)", sep=" ")
          m <- lmerTest::lmer(f1, data = rpk_norm_WLZ_frac_wlz_t)
          m_coef <- summary(m)[["coefficients"]] %>% as.data.frame() %>% rownames_to_column("coef") %>% select(-c("coef")) %>% mutate(transcript = transcript_name, linkage = linkage_name) %>% cbind(coef = c("intercept","beta1", "beta2", "beta3"))
          lmer_bulk <- lmer_bulk %>% bind_rows(m_coef)
        }
      }
      lmer_bulk <- lmer_bulk %>% rename("Pr(>|t|)"="p_value", "transcript"="tag")
      write.csv(lmer_bulk, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/rpk_norm_expressed_fract_response_to_linkage_linkageXweek_week_220819.csv")
      lmer_bulk <- read.csv("C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/rpk_norm_expressed_fract_response_to_linkage_week_linkageXweek_beta1_220818.csv")
      lmer_bulk <- lmer_bulk %>% rename("tag"="locus_tag")
      
      lmer_bulk_1 <- lmer_bulk %>% left_join(distinct(annotations))
      count(is.na(lmer_bulk_1$MAG))
      lmer_bulk_1 <- lmer_bulk_1[!is.na(lmer_bulk_1$MAG),]
      
      lmer_bulk_linkageXweek <- lmer_bulk_1 %>% filter(coef=="beta3") %>% mutate(q_value = p.adjust(p_value))
      lmer_bulk_linkage <- lmer_bulk_1 %>% filter(coef=="beta1") %>% mutate(q_value = p.adjust(p_value))
      write.csv(lmer_bulk_linkageXweek, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/rpk_norm_expressed_fract_response_to_linkage_week_linkageXweek_beta3_220822.csv")
      write.csv(lmer_bulk_linkage, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/rpk_norm_expressed_fract_response_to_linkage_week_linkageXweek_beta1_220822.csv")

#########################################################
######### model the effects of linkages on GeTMM factional expression of MAG18/19 PULs
      head(PUL_locus_tags)
      #   PUL          locus_tag        locus   tag Cons_PUL                     MAG    MAG_ID         locus_PUL
      #1 PUL1 MBINF3D6CCAE_00885 MBINF3D6CCAE 00885     <NA> P06C888_maxbin.023.mod2 MAGBg0019 MBINF3D6CCAE_PUL1
      #2 PUL1 MBINF3D6CCAE_00890 MBINF3D6CCAE 00890     <NA> P06C888_maxbin.023.mod2 MAGBg0019 MBINF3D6CCAE_PUL1
      
      #Filter GeTMM counts df to transcripts that are part of MAG18 or MAG19 named PULs
      head(rpk_norm)
      #                   P01C5551102 P01C5551105 P01C5551107 P01C6661102 P01C6661105 P01C6661107 P01C7771102 P01C7771105 P01C7771107 P01C8881102 P01C8881105 P01C8881107
      #MBINA83BCE1B_00005           0    3.944773   0.0000000           0           0           0           0           0           0           0           0           0
      #MBINA83BCE1B_00010           0    0.000000   1.7094017           0           0           0           0           0           0           0           0           0
      # Merge with MAG18/19 PUL locus tags and sum GeTMM founts by PUL
      temp <- PUL_locus_tags %>% select(c("locus_tag","locus_PUL"))
      rpk_norm_PUL_sum <- rpk_norm %>% as.data.frame %>% rownames_to_column("locus_tag") %>% inner_join(temp, by="locus_tag") %>% 
        group_by(locus_PUL) %>% summarise(across(-c("locus_tag"), sum)) %>%
        column_to_rownames("locus_PUL")
      
      rpk_norm_PUL_sum <- log2(rpk_norm_PUL_sum + 1)
      
      # transpose df, join with linkage data, center and scale
      rpk_norm_PUL_sum_t <- rpk_norm_PUL_sum %>% 
        t() %>% as.data.frame %>% rownames_to_column("SID") %>% 
        inner_join(MAG_vst_linkage_t) %>% mutate(across(-c("SID","PID","study_week"), scale))
      
      # Filter to just the 14 transcripts in fig15 and the consensus PULs
          #Make a list of consensus PULS 
          nonCons_locus_PUL <- filter(PUL_locus_tags, is.na(Cons_PUL)) %>% select("locus_PUL") %>% distinct() %>% as.list()
          Cons_locus_PUL <- filter(PUL_locus_tags, !is.na(Cons_PUL)) %>% select("locus_PUL") %>% distinct() %>% as.list()
          Cons_locus_PUL <- Cons_locus_PUL[[1]]
          Cons_locus_PUL <- Cons_locus_PUL[!Cons_locus_PUL %in% c('MBIN636B8FE1_PUL9')]
        rpk_norm_PUL_sum_t_filt <- rpk_norm_PUL_sum_t %>% select(-c(T15a_linkage_NOnames)) %>% select(-c(nonCons_locus_PUL[[1]]))
        
      # test the model: transcript VST ~ linkage quant + week + linkage quant:week + (1|PID)
      ##############################################################################
      #Where i represents sum-PUL transcripts and j represents linkages
      rm(i_length, j_length, transcript_name, linkage_name, f1, m, m_coef, lmer_bulk)
      i_length=length(Cons_locus_PUL)
      j_length=length(T15a_linkage_names)
      i=1
      j=1
      lmer_bulk <- data.frame()
      for (i in seq(1,i_length,1)) 
      {
        transcript_name=Cons_locus_PUL[i]
        for (j in seq(1,j_length,1)) 
        {
          linkage_name=T15a_linkage_names[j]
          f1 <-paste(transcript_name, "~", linkage_name, "+", "study_week", "+", linkage_name, ":study_week", "+", "(1|PID)", sep=" ")
          m <- lmerTest::lmer(f1, data = rpk_norm_PUL_sum_t_filt)
          m_coef <- summary(m)[["coefficients"]] %>% as.data.frame() %>% rownames_to_column("coef") %>% select(-c("coef")) %>% mutate(transcript = transcript_name, linkage = linkage_name) %>% cbind(coef = c("intercept","beta1", "beta2", "beta3"))
          lmer_bulk <- lmer_bulk %>% bind_rows(m_coef)
        }
      }
      lmer_bulk <- lmer_bulk %>% rename("Pr(>|t|)"="p_value", "transcript"="locus_PUL")
      write.csv(lmer_bulk, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/PUL_transcript_response_to_linkage_linkageXweek_week_220823.csv")
      #lmer_bulk <- read.csv("C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/PUL_transcript_response_to_linkage_linkageXweek_week_220823.csv") %>% select(-c("X"))
      
      PUL_locus <- PUL_locus_tags %>% select(c("locus_PUL","MAG_ID","Cons_PUL")) %>% distinct()
      lmer_bulk_1 <- lmer_bulk %>% left_join(select(PUL_locus_tags,c("locus_PUL","MAG_ID","Cons_PUL")) %>% distinct(),by="locus_PUL") %>% distinct()
      
      lmer_bulk_linkageXweek <- lmer_bulk_1 %>% filter(coef=="beta3") %>% mutate(q_value = p.adjust(p_value))
      lmer_bulk_linkage <- lmer_bulk_1 %>% filter(coef=="beta1") %>% mutate(q_value = p.adjust(p_value))
      write.csv(lmer_bulk_linkageXweek, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/PUL_transcript_response_to_linkage_linkageXweek_week_beta3_220823.csv")
      write.csv(lmer_bulk_linkage, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/PUL_transcript_response_to_linkage_linkageXweek_week_beta1__220823.csv")
      
      ##Plot
      ggplot(data  = rpk_norm_PUL_sum_t,
             aes(x = MBINFE7566FD_PUL13,
                 y = X3.4.6.Mannose,
                 color = study_week))+
        geom_point()+
        theme_minimal()
        #+ geom_line(aes(group = PID)) 
        #+ labs(title = "PUL13 expression vs. X3.4.6.Mannose quant")
      
      # test the model: transcript count in PUL ~ linkage quant + (1|PID)
      ##############################################################################
      #Where i represents sum-PUL transcripts and j represents linkages
      rm(i_length, j_length, transcript_name, linkage_name, f1, m, m_coef, lmer_bulk)
      i_length=length(Cons_locus_PUL)
      j_length=length(T15a_linkage_names)
      i=1
      j=1
      lmer_bulk <- data.frame()
      for (i in seq(1,i_length,1)) 
      {
        transcript_name=Cons_locus_PUL[i]
        for (j in seq(1,j_length,1)) 
        {
          linkage_name=T15a_linkage_names[j]
          f1 <-paste(transcript_name, "~", linkage_name, "+", "(1|PID)", sep=" ")
          m <- lmerTest::lmer(f1, data = rpk_norm_PUL_sum_t_filt)
          m_coef <- summary(m)[["coefficients"]] %>% as.data.frame() %>% rownames_to_column("coef") %>% select(-c("coef")) %>% mutate(transcript = transcript_name, linkage = linkage_name) %>% cbind(coef = c("intercept","beta1"))
          lmer_bulk <- lmer_bulk %>% bind_rows(m_coef)
        }
      }
      lmer_bulk <- lmer_bulk %>% rename("Pr(>|t|)"="p_value", "transcript"="locus_PUL")
      write.csv(lmer_bulk, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/PUL_transcript_response_to_linkage_220823.csv")
      #lmer_bulk <- read.csv("C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/PUL_transcript_response_to_linkage_linkageXweek_week_220823.csv") %>% select(-c("X"))
      
      PUL_locus <- PUL_locus_tags %>% select(c("locus_PUL","MAG_ID","Cons_PUL")) %>% distinct()
      lmer_bulk_1 <- lmer_bulk %>% left_join(select(PUL_locus_tags,c("locus_PUL","MAG_ID","Cons_PUL")) %>% distinct(),by="locus_PUL") %>% distinct()
      
      lmer_bulk_linkage <- lmer_bulk_1 %>% filter(coef=="beta1") %>% mutate(q_value = p.adjust(p_value))
      write.csv(lmer_bulk_linkage, "C:/Users/dmweb/Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/linkage_by_MAG_abnd/PUL_transcript_response_to_linkage_beta1_220823.csv")
      