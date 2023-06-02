library(tidyverse)
library(progress)
library(readxl)

datestring <- format(Sys.time(), "%y%m%d_%H%M")

setwd("~/Documents/Projects/human_studies/MDCF POC/210614_mag_annotation_database_update/")

mag_1000 <- read.table("~/Library/CloudStorage/Box-Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/MAGs/210614_0934_MDCF_POC_MAG_1000_counts_tpm_full_set.txt", header = TRUE, sep = "\t") %>%
  rownames_to_column("MAG") %>%
  select(MAG) %>%
  as.vector()

files <- list.files(path = "faa_files", pattern = "*.faa", full.names = TRUE)

bulk_table_faa <- data.frame(MAG = character(0),
                         locus_tag = character(0))

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) :eta", total = length(files))
for (i in 1:length(files)) {
  pb$tick()
  MAG <- gsub(".faa$", "", gsub(".*/NZ_", "", files[i]))
  f <- read.table(files[i], fill = TRUE, header = FALSE) %>%
    select(V1) %>%
    filter(substr(V1, 1, 1) == ">") %>%
    transmute(locus_tag = gsub(">", "", V1),
              MAG = MAG)
  bulk_table_faa <- bulk_table_faa %>%
    bind_rows(f)
}

files2 <- list.files(path = "ffn_files", pattern = "*.ffn", full.names = TRUE)

bulk_table_ffn <- data.frame(MAG = character(0),
                         locus_tag = character(0))

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) :eta", total = length(files2))
for (i in 1:length(files2)) {
  pb$tick()
  MAG <- gsub(".ffn$", "", gsub(".*/NZ_", "", files2[i]))
  f <- read.table(files2[i], fill = TRUE, header = FALSE) %>%
    select(V1) %>%
    filter(substr(V1, 1, 1) == ">") %>%
    transmute(locus_tag = gsub(">", "", V1),
              MAG = MAG)
  bulk_table_ffn <- bulk_table_ffn %>%
    bind_rows(f)
}

bulk_table <- bind_rows(bulk_table_faa, bulk_table_ffn) %>%
  mutate(MAG = gsub(".+/", "", MAG)) %>%
  distinct() %>%
  mutate(faa = ifelse(locus_tag %in% bulk_table_faa$locus_tag, "faa", NA),
         ffn = ifelse(locus_tag %in% bulk_table_ffn$locus_tag, "ffn", NA))

files_475_complete <- list.files(path = "../210524_mcSEED_aggregate_analysis_final/Annotations_2021-01/Annotation_simplified", pattern = "*", full.names = TRUE)
files_149_complete <- list.files(path = "../210524_mcSEED_aggregate_analysis_final/Annotations_2021-04/Annotation_simplified", pattern = "*", full.names = TRUE)
files_882_complete <- list.files(path = "../210524_mcSEED_aggregate_analysis_final/Annotations_2021-05/Annotation_simplified", pattern = "*", full.names = TRUE)
files_33_complete <- list.files(path = "../210614_CAZy_mcSEED_final_33/Annotation_simplified/", pattern = "*", full.names = TRUE)

files_complete <- c(files_475_complete, files_149_complete, files_882_complete, files_33_complete)

#nrow(bulk_anno_mcseed) 490391

bulk_anno_mcseed <- data.frame(MAG = character(0),
                               locus_tag = character(0))

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) :eta", total = length(files_complete))
for (i in 1:length(files_complete)) {
  pb$tick()
  MAG <- gsub("^NZ_", "", gsub(".faa$", "", gsub(".*/", "", files_complete[i])))
  f2 <- read.table(files_complete[i], sep = "\t", header = TRUE, quote="\"") %>%
    mutate(MAG = MAG,
           mcseed_strategy = "prodigal_full") %>%
    rename(locus_tag = ID) %>%
    filter(locus_tag %in% bulk_table$locus_tag)
  bulk_anno_mcseed <- bulk_anno_mcseed %>%
    bind_rows(f2)
}

files_475_partial <- list.files(path = "../210524_mcSEED_aggregate_analysis_final/Annotations_2021-01/Annotation_simplified_partial/", pattern = "*", full.names = TRUE)
files_149_partial <- list.files(path = "../210524_mcSEED_aggregate_analysis_final/Annotations_2021-04/Annotation_simplified_partial/", pattern = "*", full.names = TRUE)

files_partial <- c(files_475_partial, files_149_partial, files_882_complete)

bulk_anno_mcseed_partial <- data.frame(MAG = character(0),
                               locus_tag = character(0))

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) :eta", total = length(files_partial))
for (i in 1:length(files_partial)) {
  pb$tick()
  MAG <- gsub("^NZ_", "", gsub(".faa$", "", gsub(".*/", "", files_partial[i])))
  f3 <- read.table(files_partial[i], sep = "\t", header = TRUE, quote="\"") %>%
    mutate(MAG = MAG,
           mcseed_strategy = "prodigal_partial") %>%
    rename(locus_tag = ID) %>%
    filter(! locus_tag %in% bulk_table$locus_tag)
  if (nrow(f3) == 0) {
    next()
  }
  bulk_anno_mcseed_partial <- bulk_anno_mcseed_partial %>%
    bind_rows(f3)
}

bulk_anno_mcseed <- bulk_anno_mcseed %>%
  mutate(MAG = gsub("$", ".mod2", MAG)) %>%
  rename(mcseed_role = Role,
         mcseed_subsystem = Subsystem,
         mcseed_name = Name,
         mcseed_module1 = Module1,
         mcseed_module2 = Module2,
         mcseed_module3 = Module3)

bulk_anno_mcseed_partial <- bulk_anno_mcseed_partial %>%
  mutate(MAG = gsub("$", ".mod2", MAG)) %>%
  rename(mcseed_role = Role,
         mcseed_subsystem = Subsystem,
         mcseed_name = Name,
         mcseed_module1 = Module1,
         mcseed_module2 = Module2,
         mcseed_module3 = Module3)

sum(unique(bulk_anno_mcseed_partial$MAG) %in% unique(bulk_anno_mcseed$MAG))

bulk_anno_mcseed2 <- bind_rows(bulk_anno_mcseed, bulk_anno_mcseed_partial)

#test <- bulk_anno_mcseed2 %>%
#  filter(locus_tag %in% bulk_anno_mcseed2$locus_tag[duplicated(locus_tag)]) %>%
#  arrange(locus_tag)

#sum(duplicated(bulk_anno_mcseed2$locus_tag))

cazy_1 <- read.table("../210330_MAG_CAZyme_re-anno/MAGS_modules.txt", sep = " ", fill = TRUE, header = FALSE) %>%
  rename(cazyme_module = V2,
         locus_tag = V1) %>%
  filter(substr(locus_tag, 1, 4) == "MBIN")
cazy_2 <- read.table("../210614_CAZy_mcSEED_final_33/mat_MAGS_modules.txt", sep = " ", fill = TRUE, header = FALSE) %>%
  rename(cazyme_module = V2,
         locus_tag = V1) %>%
  filter(substr(locus_tag, 1, 4) == "MBIN")

bulk_cazy <- bind_rows(cazy_1, cazy_2) %>%
  distinct()

cazy_module_key <- read_excel("../210330_MAG_CAZyme_re-anno/MAG_results_Suzy_MCH_edit.xlsx")

bulk_cazy <- bulk_cazy %>%
  left_join(cazy_module_key, by = "cazyme_module")

bulk_table1 <- bulk_table %>%
  full_join(bulk_anno_mcseed2, by = c("MAG", "locus_tag"))

bulk_table2 <- bulk_table1 %>%
  full_join(bulk_cazy, by = "locus_tag") %>%
  filter(MAG %in% mag_1000$MAG)

sum(!is.na(bulk_table2$mcseed_subsystem))/nrow(bulk_table2)
sum(!is.na(bulk_table2$cazyme_module))/nrow(bulk_table2)

write.table(bulk_table2, paste0(datestring, "_MDCF_POC_MAGs_bulk_annotation_table.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

#bulk_table2 <- read.table("210614_1338_MDCF_POC_MAGs_bulk_annotation_table.txt", sep = "\t", header = TRUE)

bulk_table2.filt <- bulk_table2 %>%
#  filter(MAG %in% mag_1000$MAG) %>%
  filter(!is.na(mcseed_subsystem) | !is.na(cazyme_module)) %>%
#  left_join(cazy_module_key, by = "cazyme_module") %>%
  distinct()

write.table(bulk_table2.filt, paste0(datestring, "_MDCF_POC_MAGs_bulk_annotation_table_filt.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# faa ffn comp

bulk_table2 %>%
  filter(faa == "faa" & !is.na(mcseed_subsystem) | !is.na(cazyme_module)) %>%
  summarize(n()/nrow(bulk_table2))

bulk_table2 %>%
  filter(ffn == "ffn" & !is.na(mcseed_subsystem) | !is.na(cazyme_module)) %>%
  summarize(n()/nrow(bulk_table2))

bulk_table2 %>%
  filter(ffn == "ffn" & faa == "faa" & !is.na(mcseed_subsystem) | !is.na(cazyme_module)) %>%
  summarize(n()/nrow(bulk_table2))

sum(mag_1000$MAG %in% bulk_table2$MAG)

genes_from_dan <- read.table("210621_MDCF_POC_transcript_counts_raw_full_set_genes.csv", sep = "\t", quote = "") %>%
  mutate(V1 = gsub("^\"", "", V1),
         V1 = gsub("\"$", "", V1))

sum(genes_from_dan$V1 %in% bulk_table_faa$locus_tag)/nrow(genes_from_dan)
sum(genes_from_dan$V1 %in% bulk_table_ffn$locus_tag)/nrow(genes_from_dan)
sum(genes_from_dan$V1 %in% bulk_table2$locus_tag)/nrow(genes_from_dan)
sum(genes_from_dan$V1 %in% ref$locus_tag)/nrow(genes_from_dan)

###

length(unique(annotations$MAG))
length(unique(annotations.filt$MAG))
length(unique(bulk_table2$MAG))

sum(unique(annotations$MAG) %in% unique(bulk_table2$MAG))/length(unique(annotations$MAG))

length(unique(annotations.filt$locus_tag))

length(unique(bulk_table2.filt$locus_tag))

sum(unique(annotations.filt$locus_tag) %in% unique(bulk_table2.filt$locus_tag))/length(unique(annotations.filt$locus_tag))

sum(unique(annotations.filt$locus_tag) %in% unique(ref$locus_tag))/length(unique(annotations.filt$locus_tag))

write.table(annotations.filt, "test_from_dan.txt", sep = "\t")
test1 <- read.table("test_from_dan.txt", sep = "\t")
sum(unique(test1$locus_tag) %in% unique(test2$locus_tag))/length(unique(test1$locus_tag))
write.table(bulk_table2.filt, "test_from_mch.txt", sep = "\t")
test2 <- read.table("test_from_mch.txt", sep = "\t")

ref <- read.table("210706_1219_MDCF_POC_MAGs_bulk_annotation_table_filt.txt", header = TRUE, sep = "\t", quote = "\"")
sum(unique(test1$locus_tag) %in% unique(ref$locus_tag))/length(unique(test1$locus_tag))

ref2 <- read.table("~/Library/CloudStorage/Box-Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/MAGs/210706_1219_MDCF_POC_MAGs_bulk_annotation_table.txt", header = TRUE, sep = "\t", quote = "\"")

from_cyrus <- read.table("210805_cyrus_SV1_transcript_checks/CCSVD_SV1_delta_wk12_wk0_linkage_vst_transcripts.csv", sep = ",", header = TRUE)

faa <- bulk_table %>%
  filter(faa == "faa")

ffn <- bulk_table %>%
  filter(ffn == "ffn")

sum(from_cyrus$X %in% ffn$locus_tag)/nrow(from_cyrus)

missing <- read.table("210805_cyrus_SV1_transcript_checks/missingSV1transcripts.csv", sep = ",", header = TRUE) %>%
  select(X) %>%
  mutate(MAG = gsub("_.+$", "", X))

unique(missing$MAG) %in% bin_tag_key$V2

bin_tag_key <- read.table("../210608_bin_abundance_agg_and_calc_build4/210616_MDCF_POC_MAG_bin_tag_key.txt")

missing$X %in% bulk_cazy$locus_tag

bulk_cazy %>%
  filter(locus_tag %in% missing$X)

"MBIND62FFD58_06005" %in% bulk_cazy$locus_tag

"MBIND62FFD58_06005" %in% bulk_table2$locus_tag

bulk_table1 %>%
  filter(locus_tag == "MBIND62FFD58_06005")

##

nrow(ref)
sum(!ref$locus_tag %in% bulk_table2$locus_tag)

nrow(bulk_table2)
sum(!bulk_table2$locus_tag %in% ref$locus_tag)

sum(ref$locus_tag %in% bulk_table2.filt$locus_tag)/nrow(ref)
sum(bulk_table2.filt$locus_tag %in% ref$locus_tag)/nrow(bulk_table2.filt)

test <- bulk_table2.filt %>%
  filter(!locus_tag %in% ref$locus_tag)

unique(test$MAG)

test <- bulk_table2.filt %>%
  filter(!locus_tag %in% ref$locus_tag)
