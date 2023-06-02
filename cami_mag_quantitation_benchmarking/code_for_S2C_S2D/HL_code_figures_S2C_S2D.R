#Mann-Whitney test on CAMI data

library(dplyr)
library(apeglm)
library(fgsea)
library(boot)
library(RColorBrewer)
library(ggpubr)

setwd("~/Box/Gordon_Lab/Gordon_Lab/2021/CAMI_letter/genome/DASTool_Bowtie")

DAS_bowtie <- read.table("metrics_per_bin.tsv", header=TRUE, sep='\t')

colnames(DAS_bowtie) <- gsub("_bowtie", "", colnames(DAS_bowtie))

colnames(DAS_bowtie) <- gsub("__", "_", colnames(DAS_bowtie))

colnames(DAS_bowtie) <- gsub("bp_", "bp", colnames(DAS_bowtie))

setwd("~/Box/Gordon_Lab/Gordon_Lab/2021/CAMI_letter/genome/DASTool_Kallisto")

DAS_kallisto <- read.table("metrics_per_bin.tsv", header=TRUE, sep='\t')

colnames(DAS_kallisto) <- gsub("_kallisto", "", colnames(DAS_kallisto))

colnames(DAS_kallisto) <- gsub("__", "_", colnames(DAS_kallisto))

colnames(DAS_kallisto)==colnames(DAS_bowtie)

DAS_kallisto$Method <- c("DASTool_Kallisto")
DAS_bowtie$Method <- c("DASTool_Bowtie")

both_methods <- rbind(DAS_kallisto, DAS_bowtie)

colnames(both_methods)
#[1] "sample_id"                             "Bin_ID"                                "Most_abundant_genome"                 
#[4] "Purity_bp"                             "Completeness_bp"                       "Bin_size_bp"                          
#[7] "True_positives_bp"                     "True_size_of_most_abundant_genome_bp"  "Purity_seq"                           
#[10] "Completeness_seq"                      "Bin_size_seq"                          "True_positives_seq"                   
#[13] "True_size_of_most_abundant_genome_seq" "method" 

p <- ggboxplot(both_methods, x = "Method", y = "Purity_bp", color = "Method", palette = "jco", add = "jitter") + theme(legend.position = "none")
p + stat_compare_means(label = "p.signif", method = "wilcox.test", paired=FALSE) + stat_compare_means(method = "wilcox.test", paired=FALSE, label.y=0.5)

aggregate(x = both_methods$Purity_bp,                # Specify data column
          by = list(both_methods$Method),              # Specify group indicator
          FUN = mean)   

aggregate(x = both_methods$Purity_bp,                # Specify data column
          by = list(both_methods$Method),              # Specify group indicator
          FUN = sd) 

kallisto_purity <- both_methods %>%
  filter(Method=="DASTool_Kallisto") %>%
  select(Purity_bp)

p <- ggboxplot(both_methods, x = "Method", y = "Completeness_bp", color = "Method", palette = "jco", add = "jitter") + theme(legend.position = "none")
p + stat_compare_means(label = "p.signif", method = "wilcox.test", paired=FALSE) + stat_compare_means(method = "wilcox.test", paired=FALSE, label.y=1)

aggregate(x = both_methods$Completeness_bp,                # Specify data column
          by = list(both_methods$Method),              # Specify group indicator
          FUN = mean)   

aggregate(x = both_methods$Completeness_bp,                # Specify data column
          by = list(both_methods$Method),              # Specify group indicator
          FUN = sd) 

p <- ggboxplot(both_methods, x = "Method", y = "Bin_size_bp", color = "Method", palette = "jco", add = "jitter") + theme(legend.position = "none")
p + stat_compare_means(label = "p.signif", method = "wilcox.test", paired=FALSE) + stat_compare_means(method = "wilcox.test", paired=FALSE, label.y=1)

aggregate(x = both_methods$Bin_size_bp,                # Specify data column
          by = list(both_methods$Method),              # Specify group indicator
          FUN = mean)   

aggregate(x = both_methods$Bin_size_bp,                # Specify data column
          by = list(both_methods$Method),              # Specify group indicator
          FUN = sd) 


p <- ggboxplot(both_methods, x = "Method", y = "True_positives_bp", color = "Method", palette = "jco", add = "jitter") + theme(legend.position = "none")
p + stat_compare_means(label = "p.signif", method = "wilcox.test", paired=FALSE) + stat_compare_means(method = "wilcox.test", paired=FALSE, label.y=1)

aggregate(x = both_methods$True_positives_bp,                # Specify data column
          by = list(both_methods$Method),              # Specify group indicator
          FUN = mean)   

aggregate(x = both_methods$True_positives_bp,                # Specify data column
          by = list(both_methods$Method),              # Specify group indicator
          FUN = sd) 


p <- ggboxplot(both_methods, x = "Method", y = "True_size_of_most_abundant_genome_bp", color = "Method", palette = "jco", add = "jitter") + theme(legend.position = "none")
p + stat_compare_means(label = "p.signif", method = "wilcox.test", paired=FALSE) + stat_compare_means(method = "wilcox.test", paired=FALSE, label.y=1)

aggregate(x = both_methods$True_size_of_most_abundant_genome_bp,                # Specify data column
          by = list(both_methods$Method),              # Specify group indicator
          FUN = mean)   

aggregate(x = both_methods$True_size_of_most_abundant_genome_bp,                # Specify data column
          by = list(both_methods$Method),              # Specify group indicator
          FUN = sd) 


p <- ggboxplot(both_methods, x = "Method", y = "Purity_seq", color = "Method", palette = "jco", add = "jitter") + theme(legend.position = "none")
p + stat_compare_means(label = "p.signif", method = "wilcox.test", paired=FALSE) + stat_compare_means(method = "wilcox.test", paired=FALSE, label.y=1)

aggregate(x = both_methods$Purity_seq,                # Specify data column
          by = list(both_methods$Method),              # Specify group indicator
          FUN = mean)   

aggregate(x = both_methods$Purity_seq,                # Specify data column
          by = list(both_methods$Method),              # Specify group indicator
          FUN = sd) 

p <- ggboxplot(both_methods, x = "Method", y = "Completeness_seq", color = "Method", palette = "jco", add = "jitter") + theme(legend.position = "none")
p + stat_compare_means(label = "p.signif", method = "wilcox.test", paired=FALSE) + stat_compare_means(method = "wilcox.test", paired=FALSE, label.y=1)

aggregate(x = both_methods$Completeness_seq,                # Specify data column
          by = list(both_methods$Method),              # Specify group indicator
          FUN = mean)   

aggregate(x = both_methods$Completeness_seq,                # Specify data column
          by = list(both_methods$Method),              # Specify group indicator
          FUN = sd) 


p <- ggboxplot(both_methods, x = "Method", y = "Bin_size_seq", color = "Method", palette = "jco", add = "jitter") + theme(legend.position = "none")
p + stat_compare_means(label = "p.signif", method = "wilcox.test", paired=FALSE) + stat_compare_means(method = "wilcox.test", paired=FALSE, label.y=1)

aggregate(x = both_methods$Bin_size_seq,                # Specify data column
          by = list(both_methods$Method),              # Specify group indicator
          FUN = mean)   

aggregate(x = both_methods$Bin_size_seq,                # Specify data column
          by = list(both_methods$Method),              # Specify group indicator
          FUN = sd) 

p <- ggboxplot(both_methods, x = "Method", y = "True_positives_seq", color = "Method", palette = "jco", add = "jitter") + theme(legend.position = "none")
p + stat_compare_means(label = "p.signif", method = "wilcox.test", paired=FALSE) + stat_compare_means(method = "wilcox.test", paired=FALSE, label.y=1)

aggregate(x = both_methods$True_positives_seq,                # Specify data column
          by = list(both_methods$Method),              # Specify group indicator
          FUN = mean)   

aggregate(x = both_methods$True_positives_seq,                # Specify data column
          by = list(both_methods$Method),              # Specify group indicator
          FUN = sd) 

p <- ggboxplot(both_methods, x = "Method", y = "True_size_of_most_abundant_genome_seq", color = "Method", palette = "jco", add = "jitter") + theme(legend.position = "none")
p + stat_compare_means(label = "p.signif", method = "wilcox.test", paired=FALSE) + stat_compare_means(method = "wilcox.test", paired=FALSE, label.y=1)


aggregate(x = both_methods$True_size_of_most_abundant_genome_seq,                # Specify data column
          by = list(both_methods$Method),              # Specify group indicator
          FUN = mean)   

aggregate(x = both_methods$True_size_of_most_abundant_genome_seq,                # Specify data column
          by = list(both_methods$Method),              # Specify group indicator
          FUN = sd) 

complete_90_contam_5 <- both_methods %>%
  filter(Purity_bp > 0.95 & Completeness_bp > 0.90)

p <- ggboxplot(complete_90_contam_5, x = "Method", y = "Completeness_bp", color = "Method", palette = "jco", add = "jitter") + theme(legend.position = "none")
p + stat_compare_means(label = "p.signif", method = "wilcox.test", paired=FALSE) + stat_compare_means(method = "wilcox.test", paired=FALSE, label.y=1)

p <- ggboxplot(complete_90_contam_5, x = "Method", y = "Purity_bp", color = "Method", palette = "jco", add = "jitter") + theme(legend.position = "none")
p + stat_compare_means(label = "p.signif", method = "wilcox.test", paired=FALSE) + stat_compare_means(method = "wilcox.test", paired=FALSE, label.y=1)

aggregate(x = complete_90_contam_5$Completeness_bp,                # Specify data column
          by = list(complete_90_contam_5$Method),              # Specify group indicator
          FUN = mean)   

aggregate(x = complete_90_contam_5$Completeness_bp,                # Specify data column
          by = list(complete_90_contam_5$Method),              # Specify group indicator
          FUN = sd) 

aggregate(x = complete_90_contam_5$Purity_bp,                # Specify data column
          by = list(complete_90_contam_5$Method),              # Specify group indicator
          FUN = mean)   

aggregate(x = complete_90_contam_5$Purity_bp,                # Specify data column
          by = list(complete_90_contam_5$Method),              # Specify group indicator
          FUN = sd) 


complete_70_contam_5 <- both_methods %>%
  filter(Purity_bp > 0.95 & Completeness_bp > 0.70)

p <- ggboxplot(complete_70_contam_5, x = "Method", y = "Completeness_bp", color = "Method", palette = "jco", add = "jitter") + theme(legend.position = "none")
p + stat_compare_means(label = "p.signif", method = "wilcox.test", paired=FALSE) + stat_compare_means(method = "wilcox.test", paired=FALSE, label.y=1)

aggregate(x = complete_70_contam_5$Completeness_bp,                # Specify data column
          by = list(complete_70_contam_5$Method),              # Specify group indicator
          FUN = mean)   

aggregate(x = complete_70_contam_5$Completeness_bp,                # Specify data column
          by = list(complete_70_contam_5$Method),              # Specify group indicator
          FUN = sd) 


complete_50_contam_5 <- both_methods %>%
  filter(Purity_bp > 0.95 & Completeness_bp > 0.50)

aggregate(x = complete_50_contam_5$Completeness_bp,                # Specify data column
          by = list(complete_50_contam_5$Method),              # Specify group indicator
          FUN = mean)   

aggregate(x = complete_50_contam_5$Completeness_bp,                # Specify data column
          by = list(complete_50_contam_5$Method),              # Specify group indicator
          FUN = sd) 

p <- ggboxplot(complete_50_contam_5, x = "Method", y = "Completeness_bp", color = "Method", palette = "jco", add = "jitter") + theme(legend.position = "none")
p + stat_compare_means(label = "p.signif", method = "wilcox.test", paired=FALSE) + stat_compare_means(method = "wilcox.test", paired=FALSE, label.y=1)

complete_90_contam_10 <- both_methods %>%
  filter(Purity_bp > 0.90 & Completeness_bp > 0.90)

aggregate(x = complete_90_contam_10$Completeness_bp,                # Specify data column
          by = list(complete_90_contam_10$Method),              # Specify group indicator
          FUN = mean)   

aggregate(x = complete_90_contam_10$Completeness_bp,                # Specify data column
          by = list(complete_90_contam_10$Method),              # Specify group indicator
          FUN = sd) 

p <- ggboxplot(complete_90_contam_10, x = "Method", y = "Completeness_bp", color = "Method", palette = "jco", add = "jitter") + theme(legend.position = "none")
p + stat_compare_means(label = "p.signif", method = "wilcox.test", paired=FALSE) + stat_compare_means(method = "wilcox.test", paired=FALSE, label.y=1)


complete_70_contam_10 <- both_methods %>%
  filter(Purity_bp > 0.90 & Completeness_bp > 0.70)

aggregate(x = complete_70_contam_10$Completeness_bp,                # Specify data column
          by = list(complete_70_contam_10$Method),              # Specify group indicator
          FUN = mean)   

aggregate(x = complete_70_contam_10$Completeness_bp,                # Specify data column
          by = list(complete_70_contam_10$Method),              # Specify group indicator
          FUN = sd) 

p <- ggboxplot(complete_70_contam_10, x = "Method", y = "Completeness_bp", color = "Method", palette = "jco", add = "jitter") + theme(legend.position = "none")
p + stat_compare_means(label = "p.signif", method = "wilcox.test", paired=FALSE) + stat_compare_means(method = "wilcox.test", paired=FALSE, label.y=1)

complete_50_contam_10 <- both_methods %>%
  filter(Purity_bp > 0.90 & Completeness_bp > 0.50)

aggregate(x = complete_50_contam_10$Completeness_bp,                # Specify data column
          by = list(complete_50_contam_10$Method),              # Specify group indicator
          FUN = mean)   

aggregate(x = complete_50_contam_10$Completeness_bp,                # Specify data column
          by = list(complete_50_contam_10$Method),              # Specify group indicator
          FUN = sd) 

p <- ggboxplot(complete_50_contam_10, x = "Method", y = "Completeness_bp", color = "Method", palette = "jco", add = "jitter") + theme(legend.position = "none")
p + stat_compare_means(label = "p.signif", method = "wilcox.test", paired=FALSE) + stat_compare_means(method = "wilcox.test", paired=FALSE, label.y=1)



#Interesting columns:
#Purity__bp_kallisto	


#Completeness__bp_kallisto	
#Bin_size__bp_kallisto	
#True_positives__bp_kallisto	
#True_size_of_most_abundant_genome__bp_kallisto	
#Purity__seq_kallisto	
#Completeness__seq_kallisto	
#Bin_size__seq_kallisto	
#True_positives__seq_kallisto	
#True_size_of_most_abundant_genome__seq_kallisto
