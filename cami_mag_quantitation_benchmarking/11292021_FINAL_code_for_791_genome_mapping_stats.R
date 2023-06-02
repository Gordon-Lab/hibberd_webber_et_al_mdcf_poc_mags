library(viridis)
library(ggplot2)
library(dplyr)

setwd("~/2021/CAMI_dataset/791_genome_mappings/")

#Read in key connecting the anonymized and real name genomes since the kallisto counts/bowtie counts used real genome names
#and the absolute abundance files contain the anonymous naames
genome_to_id <- read.table("genome_to_id.txt")
#Change col names
colnames(genome_to_id) <- c("short_name", "long_name")


#########################################################################################
#First, going to determine false positives and negatives for kallisto samples

#Create empty data frames
false_positives_all <- c()
false_positives_all_names <- c()
mean_false_positives_all <- c()
max_false_positves_all <- c()
min_false_positives_all <- c()
false_negatives_all <- c()
pearson_correlations_all <- c()

#Start for loop through all samples
for (i in 1:length(1:64)) {
  #Make sample name
  file_name_kallisto <- paste0("sample_", i-1, "_abundance.tsv")
  #Set new working directory
  setwd("~/Box Sync/Gordon_Lab/Gordon_Lab/2021/CAMI_dataset/791_genome_mappings/all_sample_kallisto_mapping/")
  #Read in kallisto abundance files
  kallisto_abundance <- read.table(file_name_kallisto, header=TRUE)
  #Change column name to "long_name" to match with the genome_id
  colnames(kallisto_abundance)[1] <- c("long_name")
  #Remove the .fasta 
  kallisto_abundance$long_name <- gsub(".fasta", "", kallisto_abundance$long_name)
  #Merge the genome_to_id with the kallisto abundances since the real abundances are in anonymized form
  kallisto_abundance_and_names <- merge(kallisto_abundance, genome_to_id, by="long_name")
  #Remove the long name 
  kallisto_abundance_and_names$long_name <- NULL
  #Make name of sample for real abundances
  file_name_real <- paste0("abundance", i-1, ".tsv")
  #Change directory
  setwd("~/Box Sync/Gordon_Lab/Gordon_Lab/2021/CAMI_dataset/791_genome_mappings/all_sample_real_abundances/")
  #Read in the new table
  real_abundance <- read.table(file_name_real)
  #Change column names
  colnames(real_abundance) <- c("short_name", "absolute_abundance")
  #CAlculate relative abundance from absolute abundance
  real_abundance$relative_abundance_true <-  real_abundance$absolute_abundance/sum(real_abundance$absolute_abundance)
  #Merge together with the short name
  real_and_calculated_abundance <- merge(real_abundance, kallisto_abundance_and_names, by="short_name")
  #Calculate pearson correlations
  pearson_correlation <- cor(real_and_calculated_abundance$relative_abundance_true, real_and_calculated_abundance$tpm)
  #Add pearson correlation for sample to running salary
  pearson_correlations_all <- c(pearson_correlations_all, pearson_correlation)
  #Determine false positives (real value is 0, kallisto predicts >0)
  real_and_calculated_abundance$false_positive_values <- (real_and_calculated_abundance$relative_abundance_true==0 & real_and_calculated_abundance$tpm!=0)
  #isolate the genomes w false positives
  false_positives <- real_and_calculated_abundance %>%
    filter(false_positive_values==TRUE)
  #Keep track of the false positive genomes
  false_positives_all_names <- c(false_positives_all_names, false_positives$short_name)
  #Calculate the mean tpm of the false positives in the sample
  mean_false_positives <- mean(false_positives$tpm)
  #ADd in the false positives to the running tally
  mean_false_positives_all <- c(mean_false_positives_all, mean_false_positives)
  #Determine the max tpm of the false positives
  max_false_positives <- max(false_positives$tpm)
  #ADd in max to running tally
  max_false_positves_all <- c(max_false_positves_all, max_false_positives)
  #CAlculate min false positive to running tally
  min_false_positives <- min(false_positives$tpm)
  #Add min to running tally
  min_false_positives_all <- c(min_false_positives_all, min_false_positives)
  #Determine total number of false positives
  false_positives <- sum(real_and_calculated_abundance$false_positive_values)
  #Add to running tally false positives
  false_positives_all <- c(false_positives_all, false_positives)
  #Isolate false negatives (true count >0, kallisto =0)
  real_and_calculated_abundance$false_negative_values <- (real_and_calculated_abundance$relative_abundance_true!=0 & real_and_calculated_abundance$tpm==0)
  #Determine the totaal number of false negatives
  false_negatives <- sum(real_and_calculated_abundance$false_negative_values)
  #Add in false negatives to running tally
  false_negatives_all <- c(false_negatives_all, false_negatives)
}

kallisto_false_positives_all <- false_positives_all 
kallisto_false_positives_all_names  <- false_positives_all_names
kallisto_mean_false_positives_all <- mean_false_positives_all
kallisto_max_false_positves_all <- max_false_positves_all
kallisto_min_false_positives_all <- min_false_positives_all
kallisto_false_negatives_all <- false_negatives_all
kallisto_pearson_correlations_all <- pearson_correlations_all


#Make histograms of all the results
hist(pearson_correlations_all, breaks=10, main="Pearson correlation", xlab="Pearson correlation", ylab="# of samples")
hist(false_positives_all/791, breaks=10, main="% genomes false positives/sample", xlab="Number of false positives", ylab="# of samples")
hist(false_negatives_all/791, breaks=10, main="% genomes false negatives/sample", xlab="Number of false negatives", ylab="# of samples")
hist(mean_false_positives_all, breaks=10, main="Mean tpm of false positives/sample", xlab="Mean tpm of false positives/sample", ylab="# of samples")
hist(max_false_positves_all, breaks=10, main="Max tpm false positives/sample", xlab="Max tpm of false positives/sample", ylab="# of samples")
hist(min_false_positives_all, breaks=10, main="Min tpm false positives/sample", xlab="Min tpm of false positives/sample", ylab="# of samples")

##########################################################################################################

#Create a dataframe with false positive genome names
false_positives_all_names_df <- as.data.frame(false_positives_all_names)
#Determine the number of times each genome was identified as a false positive
summary_table <- as.data.frame(table(false_positives_all_names_df$false_positives_all_names))
#Order most to least
summary_table_ordered <- summary_table[order(-summary_table$Freq), ]
#Change to short name
colnames(summary_table_ordered)[1] <- c("short_name")

#Merge the summary with the genome_to_id since it has the taxa names
summary_table_ordered_w_names <- merge(summary_table_ordered, genome_to_id, by="short_name")
#Order by frequency
summary_table_ordered_w_names_ordered <- summary_table_ordered_w_names[order(-summary_table_ordered_w_names$Freq), ]
#Write this as a table
write.csv(summary_table_ordered_w_names_ordered, "frequency_of_genome_MAG_showing_up_as_false_positive_across_64_samples.csv")


setwd("~/Box Sync/Gordon_Lab/Gordon_Lab/2021/CAMI_dataset/791_genome_mappings/")
#Read in an example tax file
tax_profile <- read.table("taxonomic_profile_9.txt", header=TRUE)
#Change name to short name
colnames(tax_profile)[6] <- c("short_name")

#Make sure all 723 genomes are in the taxa profile list so that can be merged w the 723 names of the genomes
#that showed up as a false positive at least once
sum(summary_table_ordered_w_names_ordered$short_name %in% tax_profile$short_name)

#Merge the tax profile w the summary table
summary_table_ordered_w_names_tax_info <- merge(summary_table_ordered_w_names_ordered, tax_profile, by="short_name", all.x=TRUE)
#order from highest to lowest
summary_table_w_names_tax_info_ordered <- summary_table_ordered_w_names_tax_info[order(-summary_table_ordered_w_names_tax_info$Freq), ]
#remove two unimportant columns
summary_table_w_names_tax_info_ordered$PERCENTAGE <- NULL
summary_table_w_names_tax_info_ordered$TAXPATH <- NULL
#Write data to table
write.csv(summary_table_w_names_tax_info_ordered, "summary_false_positives_table_w_taxa_info_ordered.csv")

#isolate the number of times each taxonomic ID had been repeated (ie strains)
taxa_table <- as.data.frame(table(summary_table_w_names_tax_info_ordered$TAXPATHSN))
#Order
taxa_table_ordered <- taxa_table[order(-taxa_table$Freq), ]
#Write table
write.csv(taxa_table_ordered, "how_many_genomes_or_MAGs_w_same_taxa_annotation_show_up_as_false_positive.csv")

#Isolate genomes never designated as a false positive
no_false_positive <- genome_to_id %>%
  filter(!(short_name %in% summary_table_w_names_tax_info_ordered$short_name))
#Merge by the short name
no_false_positive_w_taxa <- merge(no_false_positive, tax_profile, by="short_name", all.x=TRUE)
#Remove two unimportant columns
no_false_positive_w_taxa$PERCENTAGE <- NULL
no_false_positive_w_taxa$TAXPATH <- NULL
#Write into a table
write.csv(no_false_positive_w_taxa, "table_of_taxa_with_no_false_positives.csv")

##########################################################################################################
####Now determine false positives and negatives for Bowtie samples

#Create empty data frame
false_positives_all <- c()
false_positives_all_names <- c()
mean_false_positives_all <- c()
max_false_positves_all <- c()
min_false_positives_all <- c()
false_negatives_all <- c()
pearson_correlations_all <- c()
mean_false_negatives_all <- c()
max_false_negatives_all <- c()
min_false_negatives_all <- c()

for (i in 1:length(1:64)) {
  setwd("~/Box Sync/Gordon_Lab/Gordon_Lab/2021/CAMI_dataset/791_genome_mappings/rpkm_files/")
  #Create file name
  file_name_bowtie <- paste0("sample_", i-1, "_coverage_rpkm.txt")
  #Read in table with abundances
  bowtie_abundance <- read.table(file_name_bowtie, header=TRUE, sep='\t')
  #Change column name to long name to match the real abundances later on
  colnames(bowtie_abundance)[1] <- c("long_name")
  #Remove the ending for the samples that had "fasta" rather than "fa"
  bowtie_abundance$long_name <- gsub(".fnasta", "", bowtie_abundance$long_name)
  #Remove ending mismatch in long_name
  genome_to_id$long_name <- gsub(".fasta.fasta", "", genome_to_id$long_name)
  #Remove ending mismatch in genome_to_id
  genome_to_id$long_name <- gsub(".fasta", "", genome_to_id$long_name)
  #create scaling factor to turn rpkm to tpm
  scaling_factor <- sum(as.numeric(bowtie_abundance[, 2]))/1000000
  #Convert rpkm to tpm
  bowtie_abundance[, 2] <- bowtie_abundance[, 2]/scaling_factor
  #Merge bowtie abundance with the key to compare with real abundance
  bowtie_abundance_and_names <- merge(bowtie_abundance, genome_to_id, by="long_name")
  #Remove longname
  bowtie_abundance_and_names$long_name <- NULL
  #Create name for the real abundance file
  file_name_real <- paste0("abundance", i-1, ".tsv")
  setwd("~/Box Sync/Gordon_Lab/Gordon_Lab/2021/CAMI_dataset/791_genome_mappings/all_sample_real_abundances/")
  #Read in the real abundance file
  real_abundance <- read.table(file_name_real)
  #Change the names of the columns
  colnames(real_abundance) <- c("short_name", "absolute_abundance")
  #Convert from absolute abundance to relative abundance
  real_abundance$relative_abundance_true <-  real_abundance$absolute_abundance/sum(real_abundance$absolute_abundance)
  #Merge the real and absolute abundances
  real_and_calculated_abundance <- merge(real_abundance, bowtie_abundance_and_names, by="short_name")
  #Pearson correlation between real and calculated abundances
  pearson_correlation <- cor(real_and_calculated_abundance$relative_abundance_true, real_and_calculated_abundance[, 4])
  #ADd pearson correlation to running correlations
  pearson_correlations_all <- c(pearson_correlations_all, pearson_correlation)
  #Calculate false positive (real value is zero and calculated value is non-zero)
  real_and_calculated_abundance$false_positive_values <- (real_and_calculated_abundance$relative_abundance_true==0 & real_and_calculated_abundance[, 4]!=0)
  #Filter to isolate the false positives
  false_positives <- real_and_calculated_abundance %>%
    filter(false_positive_values==TRUE)
  #Add the false positive names in
  false_positives_all_names <- c(false_positives_all_names, false_positives$short_name)
  #Calculate the mean total number of false positives per sample
  mean_false_positives <- mean(false_positives[, 4])
  #Add this to the running tally of mean false positive tpm counts in a given sample
  mean_false_positives_all <- c(mean_false_positives_all, mean_false_positives)
  #Calculate the max tpm count of all false positives in the sample
  max_false_positives <- max(false_positives[, 4])
  #Add the max false positive to the running tally
  max_false_positves_all <- c(max_false_positves_all, max_false_positives)
  #Calculate the min tpm count of all false positives in the sample
  min_false_positives <- min(false_positives[, 4])
  #Add the min false positive to the running tally
  min_false_positives_all <- c(min_false_positives_all, min_false_positives)
  #Calculate the total false positives in the sample
  false_positives <- sum(real_and_calculated_abundance$false_positive_values)
  #Add to running tally of false positives
  false_positives_all <- c(false_positives_all, false_positives)
  #Determine which genomes were false negatives
  real_and_calculated_abundance$false_negative_values <- (real_and_calculated_abundance$relative_abundance_true!=0 & real_and_calculated_abundance[, 4]==0)
  #Calculate total number of false negatives
  false_negatives <- sum(real_and_calculated_abundance$false_negative_values)
  #ADd to running tally of false negatives
  false_negatives_all <- c(false_negatives_all, false_negatives)
  if (false_negatives>0) {
    #Isolate false negatives
    false_negatives_samples <- real_and_calculated_abundance %>%
      filter(false_negative_values==TRUE)
    #Calculate mean false negative values
    mean_false_negatives <- mean(false_negatives_samples$relative_abundance_true*1000000)
    #Add to running tally mean false negatives
    mean_false_negatives_all <- c(mean_false_negatives_all, mean_false_negatives)
    #Calculate max false negative values
    max_false_negatives <- max(false_negatives_samples$relative_abundance_true*1000000)
    #Add to running tally max false negatives
    max_false_negatives_all <- c(max_false_negatives_all, max_false_negatives)
    #Calculate min false negative values
    min_false_negatives <- min(false_negatives_samples$relative_abundance_true*1000000)
    #Add to running tally min false negatives
    min_false_negatives_all <- c(min_false_negatives_all, min_false_negatives)
  }
}

bowtie_false_positives_all <- false_positives_all 
bowtie_false_positives_all_names  <- false_positives_all_names
bowtie_mean_false_positives_all <- mean_false_positives_all
bowtie_max_false_positves_all <- max_false_positves_all
bowtie_min_false_positives_all <- min_false_positives_all
bowtie_false_negatives_all <- false_negatives_all
bowtie_pearson_correlations_all <- pearson_correlations_all
bowtie_mean_false_negatives_all <- mean_false_negatives_all
bowtie_max_false_negatives_all <- max_false_negatives_all
bowtie_min_false_negatives_all <- min_false_negatives_all


false_positives_all <- as.data.frame(c(bowtie_false_positives_all, kallisto_false_positives_all))
colnames(false_positives_all) <- c("Percent_false_positives")
false_positives_all$Mapping_strategy <- c("Bowtie")
false_positives_all$Mapping_strategy[65:128] <- c("Kallisto")

write.csv(false_positives_all, "number_of_false_positives_all_samples_both_methods.csv")

ggplot(false_positives_all, aes(x=Percent_false_positives/791, fill=Mapping_strategy, color=Mapping_strategy)) + geom_histogram(alpha=0.05, bins=20, position="identity")

correlations_all <- as.data.frame(c(bowtie_pearson_correlations_all, kallisto_pearson_correlations_all))
colnames(correlations_all) <- c("Pearson_correlation")
correlations_all$Mapping_strategy <- c("Bowtie")
correlations_all$Mapping_strategy[65:128] <- c("Kallisto")

write.csv(correlations_all, "Pearson_correlations_all_64_samples.csv")

ggplot(correlations_all, aes(x=Pearson_correlation, fill=Mapping_strategy, color=Mapping_strategy)) + geom_histogram(alpha=0.1, bins=10, position="identity")

mean_false_positives_all <- as.data.frame(c(bowtie_mean_false_positives_all, kallisto_mean_false_positives_all))
colnames(mean_false_positives_all) <- c("Mean_false_positive_TPM_value")
mean_false_positives_all$Mapping_strategy <- c("Bowtie")
mean_false_positives_all$Mapping_strategy[65:128] <- c("Kallisto")

write.csv(mean_false_positives_all, "mean_tpm_count_false_positives_across_all_samples.csv")

ggplot(mean_false_positives_all, aes(x=Mean_false_positive_TPM_value, fill=Mapping_strategy, color=Mapping_strategy)) + geom_histogram(alpha=0.1, bins=20, position="identity")

max_false_positives_all <- as.data.frame(c(bowtie_max_false_positves_all, kallisto_max_false_positves_all))
colnames(max_false_positives_all) <- c("Max_false_positive_TPM_value")
max_false_positives_all$Mapping_strategy <- c("Bowtie")
max_false_positives_all$Mapping_strategy[65:128] <- c("Kallisto")

write.csv(max_false_positives_all, "max_tpm_false_positives_all_across_all_samples.csv")

ggplot(max_false_positives_all, aes(x=Max_false_positive_TPM_value, fill=Mapping_strategy, color=Mapping_strategy)) + geom_histogram(alpha=0.1, bins=20, position="identity")

min_false_positives_all <- as.data.frame(c(bowtie_min_false_positives_all, kallisto_min_false_positives_all))
colnames(min_false_positives_all) <- c("Min_false_positive_TPM_value")
min_false_positives_all$Mapping_strategy <- c("Bowtie")
min_false_positives_all$Mapping_strategy[65:128] <- c("Kallisto")

write.csv(min_false_positives_all, "min_false_positive_tpm_across_all_samples.csv")

ggplot(min_false_positives_all, aes(x=Min_false_positive_TPM_value, fill=Mapping_strategy, color=Mapping_strategy)) + geom_histogram(alpha=0.1, bins=20, position="identity")


false_negatives_all <- as.data.frame(c(bowtie_false_negatives_all,kallisto_false_negatives_all))
colnames(false_negatives_all) <- c("Percent_of_false_negatives_per_sample")
false_negatives_all$Mapping_strategy <- c("Bowtie")
false_negatives_all$Mapping_strategy[65:128] <- c("Kallisto")

write.csv(false_negatives_all, "total_number_of_false_negatives_all.csv")

ggplot(false_negatives_all, aes(x=Percent_of_false_negatives_per_sample/791, fill=Mapping_strategy, color=Mapping_strategy)) + geom_histogram(alpha=0.1, bins=20, position="identity")

bowtie_mean_false_negative_df <- as.data.frame(c(bowtie_mean_false_negatives_all))
colnames(bowtie_mean_false_negative_df) <- c("Mean_false_negative_TPM_value")
bowtie_mean_false_negative_df$Mapping_strategy <- c("Bowtie")

write.csv(bowtie_mean_false_negative_df, "mean_false_negative_tpm_value_all_samples.csv")

ggplot(bowtie_mean_false_negative_df, aes(x=Mean_false_negative_TPM_value, fill=Mapping_strategy, color=Mapping_strategy)) + geom_histogram(alpha=0.1, bins=20, position="identity")


bowtie_max_false_negative_df <- as.data.frame(c(bowtie_max_false_negatives_all))
colnames(bowtie_max_false_negative_df) <- c("Max_false_negative_TPM_value")
bowtie_max_false_negative_df$Mapping_strategy <- c("Bowtie")

write.csv(bowtie_max_false_negative_df, "max_false_negative_tpm_value_all_samples.csv")

ggplot(bowtie_max_false_negative_df, aes(x=Max_false_negative_TPM_value, fill=Mapping_strategy, color=Mapping_strategy)) + geom_histogram(alpha=0.1, bins=20, position="identity")

bowtie_min_false_negative_df <- as.data.frame(c(bowtie_min_false_negatives_all))
colnames(bowtie_min_false_negative_df) <- c("Min_false_negative_TPM_value")
bowtie_min_false_negative_df$Mapping_strategy <- c("Bowtie")

write.csv(bowtie_min_false_negative_df, "min_false_negative_tpm_values_all_samples.csv")

ggplot(bowtie_min_false_negative_df, aes(x=Min_false_negative_TPM_value, fill=Mapping_strategy, color=Mapping_strategy)) + geom_histogram(alpha=0.1, bins=20, position="identity")

#########################################################################################################
#Create histograms
hist(pearson_correlations_all, breaks=10, main="Pearson correlation", xlab="Pearson correlation", ylab="# of samples")
hist(false_positives_all/791, breaks=10, main="% genomes false positives/sample", xlab="% of false positives", ylab="# of samples")
hist(false_negatives_all/791, breaks=10, main="% genomes false negatives/sample", xlab="% of false negatives", ylab="# of samples")
hist(mean_false_positives_all, breaks=10, main="Mean tpm of false positives/sample", xlab="Mean tpm of false positives/sample", ylab="# of samples")
hist(max_false_positves_all, breaks=10, main="Max tpm false positives/sample", xlab="Max tpm of false positives/sample", ylab="# of samples")
hist(min_false_positives_all, breaks=10, main="Min tpm false positives/sample", xlab="Min tpm of false positives/sample", ylab="# of samples")
hist(mean_false_negatives_all, breaks=20, main="Mean tpm of false negatives/sample", xlab="Mean tpm of false negatives/sample", ylab="# of samples")
hist(max_false_negatives_all, breaks=20, main="Max tpm of false negatives/sample", xlab="Max tpm of false negatives/sample", ylab="# of samples")
hist(min_false_negatives_all, breaks=10, main="Min tpm of false negatives/sample", xlab="Min tpm of false negatives/sample", ylab="# of samples")

#Calculating stats

mean(bowtie_false_positives_all)
#69.3125
sd(bowtie_false_positives_all)
#28.38168
mean(kallisto_false_positives_all)
#300.1562
sd(kallisto_false_positives_all)
#50.10947
wilcox.test(bowtie_false_positives_all, kallisto_false_positives_all)
#data:  bowtie_false_positives_all and kallisto_false_positives_all
#W = 0, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0

mean(bowtie_false_negatives_all)
#17.21875
sd(bowtie_false_negatives_all)
#26.28686
mean(kallisto_false_negatives_all)
# 0.09375
sd(kallisto_false_negatives_all)
#0.4260841
wilcox.test(bowtie_false_negatives_all, kallisto_false_negatives_all)
#data:  bowtie_false_negatives_all and kallisto_false_negatives_all
#W = 4050, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0


mean(bowtie_mean_false_positives_all)
#[1] 1551.019
sd(bowtie_mean_false_positives_all)
#[1] 1760.372
mean(kallisto_mean_false_positives_all)
#[1] 22.05424
sd(kallisto_mean_false_positives_all)
#[1] 17.90934
wilcox.test(kallisto_mean_false_positives_all, bowtie_mean_false_positives_all)
#data:  kallisto_mean_false_positives_all and bowtie_mean_false_positives_all
#W = 0, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0


mean(bowtie_max_false_positves_all)
#[1] 19559.86
sd(bowtie_max_false_positves_all)
#[1] 21099.9
mean(kallisto_max_false_positves_all)
#[1] 2813.349
sd(kallisto_max_false_positves_all)
#[1] 2695.883
wilcox.test(kallisto_max_false_positves_all, bowtie_max_false_positves_all)
#data:  kallisto_max_false_positves_all and bowtie_max_false_positves_all
#W = 360, p-value = 8.844e-16
#alternative hypothesis: true location shift is not equal to 0

######################################################################################################
#Filtering analysis:

#Create empty data frames
false_positives_all <- c()
false_positives_all_names <- c()
mean_false_positives_all <- c()
max_false_positves_all <- c()
min_false_positives_all <- c()
false_negatives_all <- c()
pearson_correlations_all <- c()

for (filtered_value in 1:100) {
  for (i in 1:length(1:64)) {
    #Make sample name
    file_name_kallisto <- paste0("sample_", i-1, "_abundance.tsv")
    #Set new working directory
    setwd("~/Box Sync/Gordon_Lab/Gordon_Lab/2021/CAMI_dataset/791_genome_mappings/all_sample_kallisto_mapping/")
    #Read in kallisto abundance files
    kallisto_abundance <- read.table(file_name_kallisto, header=TRUE)
    #Change column name to "long_name" to match with the genome_id
    colnames(kallisto_abundance)[1] <- c("long_name")
    #Remove the .fasta 
    kallisto_abundance$long_name <- gsub(".fasta", "", kallisto_abundance$long_name)
    #Merge the genome_to_id with the kallisto abundances since the real abundances are in anonymized form
    kallisto_abundance_and_names <- merge(kallisto_abundance, genome_to_id, by="long_name")
    #Remove the long name 
    kallisto_abundance_and_names$long_name <- NULL
    #Make name of sample for real abundances
    file_name_real <- paste0("abundance", i-1, ".tsv")
    #Change directory
    setwd("~/Box Sync/Gordon_Lab/Gordon_Lab/2021/CAMI_dataset/791_genome_mappings/all_sample_real_abundances/")
    #Read in the new table
    real_abundance <- read.table(file_name_real)
    #Change column names
    colnames(real_abundance) <- c("short_name", "absolute_abundance")
    #CAlculate relative abundance from absolute abundance
    real_abundance$relative_abundance_true <-  real_abundance$absolute_abundance/sum(real_abundance$absolute_abundance)
    #Merge together with the short name
    real_and_calculated_abundance <- merge(real_abundance, kallisto_abundance_and_names, by="short_name")
    #Calculate pearson correlations
    pearson_correlation <- cor(real_and_calculated_abundance$relative_abundance_true, real_and_calculated_abundance$tpm)
    #Add pearson correlation for sample to running salary
    pearson_correlations_all <- c(pearson_correlations_all, pearson_correlation)
    #Determine false positives (real value is 0, kallisto predicts >0)
    real_and_calculated_abundance$false_positive_values <- (real_and_calculated_abundance$relative_abundance_true==0 & real_and_calculated_abundance$tpm>filtered_value)
    #isolate the genomes w false positives
    false_positives <- real_and_calculated_abundance %>%
      filter(false_positive_values==TRUE)
    #Keep track of the false positive genomes
    false_positives_all_names <- c(false_positives_all_names, false_positives$short_name)
    #Calculate the mean tpm of the false positives in the sample
    mean_false_positives <- mean(false_positives$tpm)
    #ADd in the false positives to the running tally
    mean_false_positives_all <- c(mean_false_positives_all, mean_false_positives)
    #Determine the max tpm of the false positives
    max_false_positives <- max(false_positives$tpm)
    #ADd in max to running tally
    max_false_positves_all <- c(max_false_positves_all, max_false_positives)
    #CAlculate min false positive to running tally
    min_false_positives <- min(false_positives$tpm)
    #Add min to running tally
    min_false_positives_all <- c(min_false_positives_all, min_false_positives)
    #Determine total number of false positives
    false_positives <- sum(real_and_calculated_abundance$false_positive_values)
    #Add to running tally false positives
    false_positives_all <- c(false_positives_all, false_positives)
    #Isolate false negatives (true count >0, kallisto =0)
    real_and_calculated_abundance$false_negative_values <- (real_and_calculated_abundance$relative_abundance_true!=0 & real_and_calculated_abundance$tpm==0)
    #Determine the totaal number of false negatives
    false_negatives <- sum(real_and_calculated_abundance$false_negative_values)
    #Add in false negatives to running tally
    false_negatives_all <- c(false_negatives_all, false_negatives)
  }
}

false_positives_all_kallisto <- false_positives_all


#Create empty data frame
false_positives_all <- c()
false_positives_all_names <- c()
mean_false_positives_all <- c()
max_false_positves_all <- c()
min_false_positives_all <- c()
false_negatives_all <- c()
pearson_correlations_all <- c()
mean_false_negatives_all <- c()
max_false_negatives_all <- c()
min_false_negatives_all <- c()

for (filtered_value in 1:100) {
  for (i in 1:length(1:64)) {
    setwd("~/Box Sync/Gordon_Lab/Gordon_Lab/2021/CAMI_dataset/791_genome_mappings/rpkm_files/")
    #Create file name
    file_name_bowtie <- paste0("sample_", i-1, "_coverage_rpkm.txt")
    #Read in table with abundances
    bowtie_abundance <- read.table(file_name_bowtie, header=TRUE, sep='\t')
    #Change column name to long name to match the real abundances later on
    colnames(bowtie_abundance)[1] <- c("long_name")
    #Remove the ending for the samples that had "fasta" rather than "fa"
    bowtie_abundance$long_name <- gsub(".fnasta", "", bowtie_abundance$long_name)
    #Remove ending mismatch in long_name
    genome_to_id$long_name <- gsub(".fasta.fasta", "", genome_to_id$long_name)
    #Remove ending mismatch in genome_to_id
    genome_to_id$long_name <- gsub(".fasta", "", genome_to_id$long_name)
    #create scaling factor to turn rpkm to tpm
    scaling_factor <- sum(as.numeric(bowtie_abundance[, 2]))/1000000
    #Convert rpkm to tpm
    bowtie_abundance[, 2] <- bowtie_abundance[, 2]/scaling_factor
    #Merge bowtie abundance with the key to compare with real abundance
    bowtie_abundance_and_names <- merge(bowtie_abundance, genome_to_id, by="long_name")
    #Remove longname
    bowtie_abundance_and_names$long_name <- NULL
    #Create name for the real abundance file
    file_name_real <- paste0("abundance", i-1, ".tsv")
    setwd("~/Box Sync/Gordon_Lab/Gordon_Lab/2021/CAMI_dataset/791_genome_mappings/all_sample_real_abundances/")
    #Read in the real abundance file
    real_abundance <- read.table(file_name_real)
    #Change the names of the columns
    colnames(real_abundance) <- c("short_name", "absolute_abundance")
    #Convert from absolute abundance to relative abundance
    real_abundance$relative_abundance_true <-  real_abundance$absolute_abundance/sum(real_abundance$absolute_abundance)
    #Merge the real and absolute abundances
    real_and_calculated_abundance <- merge(real_abundance, bowtie_abundance_and_names, by="short_name")
    #Pearson correlation between real and calculated abundances
    pearson_correlation <- cor(real_and_calculated_abundance$relative_abundance_true, real_and_calculated_abundance[, 4])
    #ADd pearson correlation to running correlations
    pearson_correlations_all <- c(pearson_correlations_all, pearson_correlation)
    #Calculate false positive (real value is zero and calculated value is non-zero)
    real_and_calculated_abundance$false_positive_values <- (real_and_calculated_abundance$relative_abundance_true==0 & real_and_calculated_abundance[, 4]>filtered_value)
    #Filter to isolate the false positives
    false_positives <- real_and_calculated_abundance %>%
      filter(false_positive_values==TRUE)
    #Add the false positive names in
    false_positives_all_names <- c(false_positives_all_names, false_positives$short_name)
    #Calculate the mean total number of false positives per sample
    mean_false_positives <- mean(false_positives[, 4])
    #Add this to the running tally of mean false positive tpm counts in a given sample
    mean_false_positives_all <- c(mean_false_positives_all, mean_false_positives)
    #Calculate the max tpm count of all false positives in the sample
    max_false_positives <- max(false_positives[, 4])
    #Add the max false positive to the running tally
    max_false_positves_all <- c(max_false_positves_all, max_false_positives)
    #Calculate the min tpm count of all false positives in the sample
    min_false_positives <- min(false_positives[, 4])
    #Add the min false positive to the running tally
    min_false_positives_all <- c(min_false_positives_all, min_false_positives)
    #Calculate the total false positives in the sample
    false_positives <- sum(real_and_calculated_abundance$false_positive_values)
    #Add to running tally of false positives
    false_positives_all <- c(false_positives_all, false_positives)
    #Determine which genomes were false negatives
    real_and_calculated_abundance$false_negative_values <- (real_and_calculated_abundance$relative_abundance_true!=0 & real_and_calculated_abundance[, 4]==0)
    #Calculate total number of false negatives
    false_negatives <- sum(real_and_calculated_abundance$false_negative_values)
    #ADd to running tally of false negatives
    false_negatives_all <- c(false_negatives_all, false_negatives)
    if (false_negatives>0) {
      #Isolate false negatives
      false_negatives_samples <- real_and_calculated_abundance %>%
        filter(false_negative_values==TRUE)
      #Calculate mean false negative values
      mean_false_negatives <- mean(false_negatives_samples$relative_abundance_true*1000000)
      #Add to running tally mean false negatives
      mean_false_negatives_all <- c(mean_false_negatives_all, mean_false_negatives)
      #Calculate max false negative values
      max_false_negatives <- max(false_negatives_samples$relative_abundance_true*1000000)
      #Add to running tally max false negatives
      max_false_negatives_all <- c(max_false_negatives_all, max_false_negatives)
      #Calculate min false negative values
      min_false_negatives <- min(false_negatives_samples$relative_abundance_true*1000000)
      #Add to running tally min false negatives
      min_false_negatives_all <- c(min_false_negatives_all, min_false_negatives)
    }
  }
}

false_positives_all_bowtie <- false_positives_all

#Filter 1
wilcox.test(false_positives_all_kallisto[1:64], false_positives_all_bowtie[1:64])
#p-value = 7.344e-13

#Filter 2
wilcox.test(false_positives_all_kallisto[65:128], false_positives_all_bowtie[65:128])
#p-value = 4.598e-05

#Filter 3
wilcox.test(false_positives_all_kallisto[129:192], false_positives_all_bowtie[129:192])
#p-value = p-value = 0.344

#Filter 4
wilcox.test(false_positives_all_kallisto[193:256], false_positives_all_bowtie[193:256])
#p-value = p-value = 0.2617

#Filter 5
wilcox.test(false_positives_all_kallisto[257:320], false_positives_all_bowtie[257:320])
#p-value = p-value = 0.02013

