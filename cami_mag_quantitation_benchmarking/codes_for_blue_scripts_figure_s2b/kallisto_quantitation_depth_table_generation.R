####Code for kallisto_quantification_depth_table_generation.R######

library(progress)
library(rhdf5)
library(matrixStats)
library(dplyr)
library(tibble)

setwd("/HL_working_directory/03292021_abundances_for_anonymous_contigs/all_samples_h5_from_kallisto")

tpm <- function(counts, effLen) {
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1E6))
}

#subgroup_file <- read.table("subgroup_file.txt", sep = "\t")
#colnames(subgroup_file) <- c("SID", "PID")

#participants <- c("sample_0", "sample_1", "sample_2", "sample_3", "sample_4", "sample_5", "sample_6", "sample_7", "sample_8", "sample_9")

#pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) :eta", total = 952)
#for (i in 1:length(participants)) {

files <- list.files(".", "*.h5", full.names = TRUE, recursive = TRUE)

SID_array <- NA

res <- data.frame()

for (j in 1:length(files)) {
  #    pb$tick()
  SID <- gsub(".h5$", "", basename(files[j]))
  if (is.na(SID_array[1])) {
    SID_array <- SID
  } else {
    SID_array <- c(SID_array, SID)
  }
  boot_array <- h5read(files[j], "bootstrap")
  df_boot <- data.frame(contigName = h5read(files[j], "aux/ids"),
                        contigLen = h5read(files[j], "/aux/lengths"),
                        contigEffLen = h5read(files[j], "/aux/eff_lengths"))
  
  if (length(res) == 0) {
    res <- df_boot %>%
      select(-contigEffLen)
  }
  
  for (k in 1:length(boot_array)) {
    df_boot <- df_boot %>%
      bind_cols(data.frame(tpm(boot_array[[k]], df_boot$contigEffLen)))
    #%>% rename_with(function(x) paste0("bs", k - 1)))
    colnames(df_boot)[k+3] <- paste0("bs", k-1)
  }
  
  res_kallisto <- read.table(list.files(".", paste0(SID, ".tsv"), full.names = TRUE, recursive = TRUE), sep = "\t", header = TRUE)
  
  if (!all(df_boot$contigName == res_kallisto$target_id)) {
    stop()
  }
  
  res_sub <- data.frame(ab = res_kallisto$tpm,
                        #                          mean_boot = rowMeans(df_boot %>% select(-contigName, -contigLen, -contigEffLen)),
                        var = rowVars(as.matrix(df_boot %>% select(-contigName, -contigLen, -contigEffLen))),
                        contigName = df_boot$contigName,
                        contigLen = df_boot$contigLen)
  colnames(res_sub)[colnames(res_sub) == "ab"] <- SID
  colnames(res_sub)[colnames(res_sub) == "var"] <- paste0(SID, "-var")
  
  res <- res %>%
    full_join(res_sub, by = c("contigName", "contigLen"))
  
}

# res <- res %>%
#   select(SID_array) %>%
#   transmute(totalAvgDepth = rowSums(.)) %>%
#   bind_cols(res) %>%
#   relocate(totalAvgDepth, .after = contigLen)

all_depths <- res %>%
  select(all_of(SID_array))
totalAvgDepth <- rowSums(all_depths)
res_final <- add_column(res, totalAvgDepth, .after="contigLen")

write.table(res_final, paste0("CAMI_mouse_to_anonymous_contigs_kallisto.depth"), sep = "\t", quote = F, row.names = F)