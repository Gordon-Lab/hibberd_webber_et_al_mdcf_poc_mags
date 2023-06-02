library(viridis)
library(ggplot2)
library(dplyr)

setwd("~/Box Sync/Gordon_Lab/Gordon_Lab/2021/CAMI_dataset/791_genome_mappings/")

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


#kallisto_false_positives_all <- false_positives_all 
#kallisto_false_positives_all_names  <- false_positives_all_names
#kallisto_mean_false_positives_all <- mean_false_positives_all
#kallisto_max_false_positves_all <- max_false_positves_all
#kallisto_min_false_positives_all <- min_false_positives_all
#kallisto_false_negatives_all <- false_negatives_all
#kallisto_pearson_correlations_all <- pearson_correlations_all

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

#Calculating stats for POC paper 

mean(bowtie_mean_false_positives_all)
#[1] 1551.019
mean(kallisto_mean_false_positives_all)
#[1] 22.05424
sd(bowtie_mean_false_positives_all)
#[1] 1760.372
sd(kallisto_mean_false_positives_all)
#[1] 17.90934
wilcox.test(kallisto_mean_false_positives_all, bowtie_mean_false_positives_all)
#p-value < 2.2e-16

#Wilcoxon rank sum test with continuity correction

#data:  kallisto_max_false_positves_all and bowtie_max_false_positves_all
#W = 360, p-value = 8.844e-16
#alternative hypothesis: true location shift is not equal to 0

mean(bowtie_max_false_positves_all)
#[1] 19559.86
mean(kallisto_max_false_positves_all)
#[1] 2813.349
sd(bowtie_max_false_positves_all)
#[1] 21099.9
sd(kallisto_max_false_positves_all)
#[1] 2695.883
wilcox.test(kallisto_max_false_positves_all, bowtie_max_false_positves_all)

#data:  kallisto_max_false_positves_all and bowtie_max_false_positves_all
#W = 360, p-value = 8.844e-16
#alternative hypothesis: true location shift is not equal to 0

##############################################################################################

#Create a table with all 64 sample abundances from kallisto so I can determine false negatives/positives at different
#tpm cutoffs

#Isolate the names of the genomes
kallisto_tpm_table <- as.data.frame(genome_to_id$short_name)
colnames(kallisto_tpm_table) <- c("short_name")


for (i in 1:length(1:64)) {
  #Create the sample names
  file_name_kallisto <- paste0("sample_", i-1, "_abundance.tsv")
  #change working directory
  setwd("~/Box Sync/Gordon_Lab/Gordon_Lab/2021/CAMI_dataset/791_genome_mappings/all_sample_kallisto_mapping/")
  #REad in kallisto abundance file
  kallisto_abundance <- read.table(file_name_kallisto, header=TRUE)
  #Change column name
  colnames(kallisto_abundance)[1] <- c("long_name")
  #REmove any remaining .fasta
  kallisto_abundance$long_name <- gsub(".fasta", "", kallisto_abundance$long_name)
  #Isolate two columns
  kallisto_tpm <- kallisto_abundance %>%
    select(c("long_name", "tpm"))
  #Re-name the tpm column w the sample name
  colnames(kallisto_tpm)[2] <- paste0("sample_", i-1, "_tpm")
  #Merge together with the genome names in order to get short name
  kallisto_tpm_and_names <- merge(kallisto_tpm, genome_to_id, by="long_name")
  #Remove long name
  kallisto_tpm_and_names$long_name <- NULL
  #Merge the individual sample counts into the larger table
  kallisto_tpm_table <- merge(kallisto_tpm_table, kallisto_tpm_and_names, by="short_name")
}

#Name the rows of the dataframe
rownames(kallisto_tpm_table) <- kallisto_tpm_table$short_name
#Remove the short_name column so all that is left are the counts
kallisto_tpm_table$short_name <- NULL
#Convert to a matrix
kallisto_tpm_table <- as.matrix(kallisto_tpm_table)

##################################################################################################

#Now going to do the same thing for the real abundances
#Create a data frame containing the short names of the genomes
real_abundance_table <- as.data.frame(genome_to_id$short_name)
colnames(real_abundance_table) <- c("short_name")

for (i in 1:length(1:64)) {
  #Create name of sample file
  file_name_real <- paste0("abundance", i-1, ".tsv")
  #Change working directory
  setwd("~/Box Sync/Gordon_Lab/Gordon_Lab/2021/CAMI_dataset/791_genome_mappings/all_sample_real_abundances/")
  #Read real abundance table
  real_abundance <- read.table(file_name_real)
  #Change column names
  colnames(real_abundance) <- c("short_name", "absolute_abundance")
  #Convert absolute abundance into relative abundance * 1000000 to mimic tpm
  real_abundance$relative_abundance_true <-  (real_abundance$absolute_abundance/sum(real_abundance$absolute_abundance))*1000000
  #REmove absolute abundance file
  real_abundance$absolute_abundance <- NULL
  #Add in sample name to column
  colnames(real_abundance)[2] <- paste0("sample_", i-1, "_relative_abundance_true")
  #Add column to exisiting table
  real_abundance_table <- merge(real_abundance_table, real_abundance, by="short_name")
}

#Add in row names
rownames(real_abundance_table) <-  real_abundance_table$short_name
#Remove the short name column
real_abundance_table$short_name <- NULL
#convert to a matrix
real_abundance_table <- as.matrix(real_abundance_table)

###########################################################################################################

#Now going to perform a sweep to ID the rates of false positives and negatives at different filtering criteria

#Create a dataframe containing important info
kallisto_abundance_tpm_filt.stat <- data.frame(threshold = numeric(0), 
                                         prevalence = numeric(0),
                                         ntaxa = numeric(0),
                                         fraction = numeric(0),
                                         false_negatives=numeric(0),
                                         false_positives=numeric(0),
                                         pos_to_neg_ratio=numeric(0))

for (threshold in seq(0, 100, 5)) {
  for (prevalence in seq(0.1, 1, 0.1)) {
    #GEnomes to keep are genomes above a certain tpm count in a give percentage of samples
    genomes_to_keep <- rowSums(kallisto_tpm_table > threshold) > (ncol(kallisto_tpm_table) * prevalence)
    #Genomes to be set to zero are the one that do not pass that filtering criteria (ie counts would in theory be removed)
    genomes_to_zero <- (rowSums(kallisto_tpm_table > threshold) <= (ncol(kallisto_tpm_table) * prevalence))
    #Creating a copy of the table to set the genomes not passing filtering to zero
    modified_kallisto_tpm_table <- kallisto_tpm_table
    #Set filtered genome counts to zero
    modified_kallisto_tpm_table[genomes_to_zero, ] <- 0
    #Identify the false positives at the filtering criteria (ie modified kallisto will be non-zero, real abundance=0)
    false_positives <- (modified_kallisto_tpm_table > 0 & real_abundance_table==0)
    #Determine the total number of false positives across all samples (out of 791 genomes*64 samples)
    false_positives_all <- sum(false_positives)
    #Identify false negatives at filtering criteria (ie modified kallisto will be zero when real abundance=0)
    false_negatives <- (modified_kallisto_tpm_table==0 & real_abundance_table>0)
    #Determine total number of false negatives across all samples (out of 791 genomes *64 samples)
    false_negatives_all <- sum(false_negatives)
    #Filtered table is the one with filtered genome rows removed
    kallisto_abundance_tpm.filt <- kallisto_tpm_table[rowSums(kallisto_tpm_table > threshold) > (ncol(kallisto_tpm_table) * prevalence),]
    #update the filtered table stats:
    kallisto_abundance_tpm_filt.stat <- bind_rows(kallisto_abundance_tpm_filt.stat,
                                            data.frame(threshold = threshold,
                                                       prevalence = prevalence,
                                                       false_positives=(false_positives_all/50624)*100,
                                                       false_negatives=(false_negatives_all/50624)*100,
                                                       pos_to_neg_ratio=(false_positives_all)/(false_negatives_all),
                                                       ntaxa = nrow(kallisto_abundance_tpm.filt),
                                                       fraction = ifelse(nrow(kallisto_abundance_tpm.filt) > 0, sum(kallisto_abundance_tpm.filt)/sum(kallisto_tpm_table), 0)))
  }
}

write.csv(kallisto_abundance_tpm_filt.stat, "kallisto_abundance_tpm_filt_stat.csv")

###############################################################################################################
#Now, make plots w the results:

#Making the ratio log so that the color looks better:
kallisto_abundance_tpm_filt.stat_w_log <- kallisto_abundance_tpm_filt.stat
kallisto_abundance_tpm_filt.stat_w_log$log_pos_to_neg_ratio <- log(kallisto_abundance_tpm_filt.stat_w_log$pos_to_neg_ratio + 1)


#Making plots
ggplot(kallisto_abundance_tpm_filt.stat_w_log, aes(x = prevalence, y = threshold, label = round((log_pos_to_neg_ratio), 2))) +
  geom_raster(aes(fill = log_pos_to_neg_ratio)) +
  scale_x_continuous(breaks = seq(0.1, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  theme_bw() +
  scale_fill_gradient2(low="yellow", high="blue", guide="colorbar")+
  ggtitle("Ratio of false positives to false negatives") 

ggplot(kallisto_abundance_tpm_filt.stat, aes(x = prevalence, y = threshold, label = round(pos_to_neg_ratio, 2))) +
  ggtitle("Ratio of false positives to false negatives") +
  geom_tile(color = "black", fill = NA) +
  geom_text() +
  scale_x_continuous(breaks = seq(0.1, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  theme_bw()

ggplot(kallisto_abundance_tpm_filt.stat, aes(x = prevalence, y = threshold, label = round((false_negatives), 2))) +
  geom_raster(aes(fill = false_negatives)) +
  scale_x_continuous(breaks = seq(0.1, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  theme_bw() +
  scale_fill_gradient2(low="yellow", high="blue", guide="colorbar")+
  ggtitle("Percent of false negatives remaining") 

ggplot(kallisto_abundance_tpm_filt.stat, aes(x = prevalence, y = threshold, label = round(false_negatives, 2))) +
  ggtitle("Percent of false negatives remaining") +
  geom_tile(color = "black", fill = NA) +
  geom_text() +
  scale_x_continuous(breaks = seq(0.1, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  theme_bw()


ggplot(kallisto_abundance_tpm_filt.stat, aes(x = prevalence, y = threshold, label = round((false_positives), 2))) +
  geom_raster(aes(fill = false_positives)) +
  scale_x_continuous(breaks = seq(0.1, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  theme_bw() +
  scale_fill_gradient2(low="yellow", high="blue", guide="colorbar")+
  ggtitle("Percent of false positives remaining") 

ggplot(kallisto_abundance_tpm_filt.stat, aes(x = prevalence, y = threshold, label = round(false_positives, 2))) +
  ggtitle("Percent of false positives remaining") +
  geom_tile(color = "black", fill = NA) +
  geom_text() +
  scale_x_continuous(breaks = seq(0.1, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  theme_bw()

write.csv(kallisto_abundance_tpm_filt.stat, 'parameter_sweep_kallisto_abundance_tpm_filt_stat.csv')

ggplot(kallisto_abundance_tpm_filt.stat, aes(x = prevalence, y = threshold, label = round((ntaxa/791), 2))) +
  geom_raster(aes(fill = ntaxa)) +
  scale_x_continuous(breaks = seq(0.1, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  theme_bw() +
  scale_fill_gradient2(low="yellow", high="blue", guide="colorbar")+
  ggtitle("Total taxa remaining") 

ggplot(kallisto_abundance_tpm_filt.stat, aes(x = prevalence, y = threshold, label = round(ntaxa/791, 2))) +
  ggtitle("Number of taxa remaining") +
  geom_tile(color = "black", fill = NA) +
  geom_text() +
  scale_x_continuous(breaks = seq(0.1, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  theme_bw()

ggplot(kallisto_abundance_tpm_filt.stat, aes(x = prevalence, y = threshold, label = round((ntaxa), 2))) +
  geom_raster(aes(fill = ntaxa)) +
  scale_x_continuous(breaks = seq(0.1, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  theme_bw() +
  scale_fill_gradient2(low="yellow", high="blue", guide="colorbar")+
  ggtitle("Total taxa remaining")

ggplot(kallisto_abundance_tpm_filt.stat, aes(x = prevalence, y = threshold, label = round((fraction), 2))) +
  geom_raster(aes(fill = fraction)) +
  scale_x_continuous(breaks = seq(0.1, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  theme_bw() +
  scale_fill_gradient2(low="yellow", high="blue", guide="colorbar")+
  ggtitle("Percent of reads remaining")


ggplot(kallisto_abundance_tpm_filt.stat, aes(x = prevalence, y = threshold, label = round(fraction, 2))) +
  ggtitle("Fraction of total counts remaining") +
  geom_tile(color = "black", fill = NA) +
  geom_text() +
  scale_x_continuous(breaks = seq(0.1, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  theme_bw()

##########################################################################################################

MAGs_to_keep <- rowSums(kallisto_tpm_table>20) > (64*0.3)
filtered_kallisto_tpm_table <- kallisto_tpm_table[MAGs_to_keep, ]
#236 MAGs by 64 samples

a <- as.data.frame(rownames(filtered_kallisto_tpm_table))
colnames(a) <- c("short_name")

b <- merge(a, genome_to_id, by="short_name", all.x=TRUE)

genomes_to_keep_onward <- b$long_name

write.csv(genomes_to_keep_onward, "names_to_keep.csv")


#####################################################################################

#Now, repeat for the 236 mapping

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
  setwd("~/Box Sync/Gordon_Lab/Gordon_Lab/2021/CAMI_dataset/236_all_sample_kallisto_mapping/")
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

kallisto_236_false_positives_all <- false_positives_all 
kallisto_236_false_positives_all_names  <- false_positives_all_names
kallisto_236_mean_false_positives_all <- mean_false_positives_all
kallisto_236_max_false_positves_all <- max_false_positves_all
kallisto_236_min_false_positives_all <- min_false_positives_all
kallisto_236_false_negatives_all <- false_negatives_all
kallisto_236_pearson_correlations_all <- pearson_correlations_all


#Make histograms of all the results
hist(kallisto_236_pearson_correlations_all, breaks=10, main="Pearson correlation", xlab="Pearson correlation", ylab="# of samples")
hist(kallisto_236_false_positives_all/236, breaks=10, main="% genomes false positives/sample", xlab="Number of false positives", ylab="# of samples")
hist(kallisto_236_false_negatives_all/236, breaks=1, main="% genomes false negatives/sample", xlab="Number of false negatives", ylab="# of samples")
hist(kallisto_236_mean_false_positives_all, breaks=10, main="Mean tpm of false positives/sample", xlab="Mean tpm of false positives/sample", ylab="# of samples")
hist(kallisto_236_max_false_positves_all, breaks=10, main="Max tpm false positives/sample", xlab="Max tpm of false positives/sample", ylab="# of samples")
hist(kallisto_236_min_false_positives_all, breaks=10, main="Min tpm false positives/sample", xlab="Min tpm of false positives/sample", ylab="# of samples")

##########################################################################################

false_positives_all <- as.data.frame(c(bowtie_false_positives_all, kallisto_false_positives_all, kallisto_236_false_positives_all))
colnames(false_positives_all) <- c("Percent_false_positives")
false_positives_all$Mapping_strategy <- c("Bowtie")
false_positives_all$Mapping_strategy[65:128] <- c("Kallisto_791")
false_positives_all$Mapping_strategy[129:192] <- c("Kallisto_236")
false_positives_all$Percent_false_positives[1:128] <- false_positives_all$Percent_false_positives[1:128]/791
false_positives_all$Percent_false_positives[129:192] <- false_positives_all$Percent_false_positives[129:192]/236

ggplot(false_positives_all, aes(x=Percent_false_positives, fill=Mapping_strategy, color=Mapping_strategy)) + geom_histogram(alpha=0.05, bins=20, position="identity")

correlations_all <- as.data.frame(c(bowtie_pearson_correlations_all, kallisto_pearson_correlations_all, kallisto_236_pearson_correlations_all))
colnames(correlations_all) <- c("Pearson_correlation")
correlations_all$Mapping_strategy <- c("Bowtie")
correlations_all$Mapping_strategy[65:128] <- c("Kallisto_791")
correlations_all$Mapping_strategy[129:192] <- c("Kallisto_236")

ggplot(correlations_all, aes(x=Pearson_correlation, fill=Mapping_strategy, color=Mapping_strategy)) + geom_histogram(alpha=0.1, bins=10, position="identity")

mean_false_positives_all <- as.data.frame(c(bowtie_mean_false_positives_all, kallisto_mean_false_positives_all, kallisto_236_mean_false_positives_all))
colnames(mean_false_positives_all) <- c("Mean_false_positive_TPM_value")
mean_false_positives_all$Mapping_strategy <- c("Bowtie")
mean_false_positives_all$Mapping_strategy[65:128] <- c("Kallisto_791")
mean_false_positives_all$Mapping_strategy[129:192] <- c("Kallisto_236")

ggplot(mean_false_positives_all, aes(x=Mean_false_positive_TPM_value, fill=Mapping_strategy, color=Mapping_strategy)) + geom_histogram(alpha=0.1, bins=20, position="identity")

max_false_positives_all <- as.data.frame(c(bowtie_max_false_positves_all, kallisto_max_false_positves_all, kallisto_236_max_false_positves_all))
colnames(max_false_positives_all) <- c("Max_false_positive_TPM_value")
max_false_positives_all$Mapping_strategy <- c("Bowtie")
max_false_positives_all$Mapping_strategy[65:128] <- c("Kallisto_791")
max_false_positives_all$Mapping_strategy[129:192] <- c("Kallisto_236")

ggplot(max_false_positives_all, aes(x=Max_false_positive_TPM_value, fill=Mapping_strategy, color=Mapping_strategy)) + geom_histogram(alpha=0.1, bins=20, position="identity")

min_false_positives_all <- as.data.frame(c(bowtie_min_false_positives_all, kallisto_min_false_positives_all, kallisto_236_min_false_positives_all))
colnames(min_false_positives_all) <- c("Min_false_positive_TPM_value")
min_false_positives_all$Mapping_strategy <- c("Bowtie")
min_false_positives_all$Mapping_strategy[65:128] <- c("Kallisto_791")
min_false_positives_all$Mapping_strategy[129:192] <- c("Kallisto_236")

ggplot(min_false_positives_all, aes(x=Min_false_positive_TPM_value, fill=Mapping_strategy, color=Mapping_strategy)) + geom_histogram(alpha=0.1, bins=20, position="identity")

false_negatives_all <- as.data.frame(c(bowtie_false_negatives_all,kallisto_false_negatives_all, kallisto_236_false_negatives_all))
colnames(false_negatives_all) <- c("Percent_of_false_negatives_per_sample")
false_negatives_all$Mapping_strategy <- c("Bowtie")
false_negatives_all$Mapping_strategy[65:128] <- c("Kallisto_791")
false_negatives_all$Mapping_strategy[129:192] <- c("Kallisto_236")
false_negatives_all$Percent_of_false_negatives_per_sample[1:128] <- false_negatives_all$Percent_of_false_negatives_per_sample[1:128]/791
false_negatives_all$Percent_of_false_negatives_per_sample[129:192] <- false_negatives_all$Percent_of_false_negatives_per_sample[129:192]/236

ggplot(false_negatives_all, aes(x=Percent_of_false_negatives_per_sample/791, fill=Mapping_strategy, color=Mapping_strategy)) + geom_histogram(alpha=0.1, bins=20, position="identity")


bowtie_mean_false_negative_df <- as.data.frame(c(bowtie_mean_false_negatives_all))
colnames(bowtie_mean_false_negative_df) <- c("Mean_false_negative_TPM_value")
bowtie_mean_false_negative_df$Mapping_strategy <- c("Bowtie")

ggplot(bowtie_mean_false_negative_df, aes(x=Mean_false_negative_TPM_value, fill=Mapping_strategy, color=Mapping_strategy)) + geom_histogram(alpha=0.1, bins=20, position="identity")


bowtie_max_false_negative_df <- as.data.frame(c(bowtie_max_false_negatives_all))
colnames(bowtie_max_false_negative_df) <- c("Max_false_negative_TPM_value")
bowtie_max_false_negative_df$Mapping_strategy <- c("Bowtie")

ggplot(bowtie_max_false_negative_df, aes(x=Max_false_negative_TPM_value, fill=Mapping_strategy, color=Mapping_strategy)) + geom_histogram(alpha=0.1, bins=20, position="identity")

bowtie_min_false_negative_df <- as.data.frame(c(bowtie_min_false_negatives_all))
colnames(bowtie_min_false_negative_df) <- c("Min_false_negative_TPM_value")
bowtie_min_false_negative_df$Mapping_strategy <- c("Bowtie")

ggplot(bowtie_min_false_negative_df, aes(x=Min_false_negative_TPM_value, fill=Mapping_strategy, color=Mapping_strategy)) + geom_histogram(alpha=0.1, bins=20, position="identity")


ggplot(kallisto_abundance_tpm_filt.stat, aes(x = prevalence, y = threshold, label = round(false_negatives, 2))) +
  ggtitle("False negatives remaining") +
  geom_tile(color = "black", fill = NA) +
  geom_text() +
  scale_x_continuous(breaks = seq(0.1, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  theme_bw()

ggplot(kallisto_abundance_tpm_filt.stat, aes(x = prevalence, y = threshold, label = round(pos_to_neg_ratio, 2))) +
  ggtitle("False positive to false negative ratio") +
  geom_tile(color = "black", fill = NA) +
  geom_text() +
  scale_x_continuous(breaks = seq(0.1, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  theme_bw()

ggplot(kallisto_abundance_tpm_filt.stat, aes(x = prevalence, y = threshold, label = round(ntaxa, 2))) +
  ggtitle("False negatives remaining") +
  geom_tile(color = "black", fill = NA) +
  geom_text() +
  scale_x_continuous(breaks = seq(0.1, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  theme_bw()

ggplot(kallisto_abundance_tpm_filt.stat, aes(x = prevalence, y = threshold, label = round(fraction, 2))) +
  ggtitle("False negatives remaining") +
  geom_tile(color = "black", fill = NA) +
  geom_text() +
  scale_x_continuous(breaks = seq(0.1, 1, 0.1)) +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  theme_bw()


x <- c(real_abundance_table$sample_0_relative_abundance_true==0 & kallisto_tpm_table$sample_0_tpm!=0)


#going to run through kallisto_abundance
#sub das_kallisto_tpm for kallisto_abundance

as.numeric(kallisto_abundance)

kallisto_abundance_filt.stat <- data.frame(threshold = numeric(0), 
                                         prevalence = numeric(0),
                                         ntaxa = numeric(0),
                                         fraction = numeric(0))

for (threshold in seq(0, 100, 5)) {
  for (prevalence in seq(0.1, 1, 0.1)) {
    kallisto_abundance.filt <- kallisto_abundance[rowSums(kallisto_abundance > threshold) > (ncol(kallisto_abundance) * prevalence),]
    kallisto_abundance_filt.stat <- bind_rows(kallisto_abundance_filt.stat,
                                            data.frame(threshold = threshold,
                                                       prevalence = prevalence,
                                                       ntaxa = nrow(kallisto_abundance.filt),
                                                       fraction = ifelse(nrow(kallisto_abundance.filt) > 0, sum(kallisto_abundance.filt)/sum(kallisto_abundance), 0)))
  }
}


