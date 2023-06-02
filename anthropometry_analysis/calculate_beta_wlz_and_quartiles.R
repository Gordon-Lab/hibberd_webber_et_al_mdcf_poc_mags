library(readxl)

## calculate beta-WLZ for each child

setwd("")

datestring <- format(Sys.time(), "%y%m%d_%H%M")

map <- read_excel("mapping_file.xlsx") %>%
  filter(sample_type == "fecal")

map.anthropometry <- read.table("MDCF_POC_anthropometry.txt", sep = "\t", header = T)
map.anthropometry$PID <- factor(map.anthropometry$PID)

map_merge <- data.frame(SID = colnames(dat_kallisto_vst)) %>%
  merge(map %>% select(SID, PID, study_week, study_phase_alt), by = "SID") %>%
  merge(map.anthropometry %>% select(-SID), by = c("PID", "study_week")) %>%
  arrange(SID)

beta_table <- data.frame(PID = character(0), coef = numeric(0), delta_wlz = numeric(0), study_arm = character(0))

for (i in 1:length(levels(map.anthropometry$PID))) {
  dat_sub <- map.anthropometry %>%
    filter(PID == levels(map.anthropometry$PID)[i] & study_phase != "wash_out") %>%
    arrange(age_at_sampling)
  
  if (any(is.na(dat_sub$wlz))) {
    stop(paste0("[STATUS] Missing or incomplete data for sample: ", levels(map.anthropometry$PID)[i]))
  } else if (nrow(dat_sub) < 3) {
    next()
  }
  
  model <- lm(wlz ~ as.numeric(as.character(study_week)), data = dat_sub)
  #  model <- lm(wlz ~ age_at_sampling, data = dat_sub)
  
  beta_table <- rbind(beta_table, data.frame(PID = levels(map.anthropometry$PID)[i], 
                                             study_arm = unique(dat_sub$study_arm),
                                             coef = unname(model$coefficients[2]),
                                             delta_wlz = dat_sub$wlz[length(dat_sub$wlz)] - dat_sub$wlz[1]))
  
}

#beta_table$beta_quartile <- cut(beta_table$coef, breaks = quantile(beta_table$coef), labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest=TRUE)
#beta_table <- beta_table %>%
#  mutate(beta_half = ifelse(beta_quartile %in% c("Q3", "Q4"), "Q3Q4", "Q1Q2"))
#beta_table$delta_quartile <- cut(beta_table$delta_wlz, breaks = quantile(beta_table$delta_wlz), labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest=TRUE)

beta_table_v2 <- beta_table %>%
  group_by(study_arm) %>%
  mutate(beta_quartile = cut(coef, breaks = quantile(coef), labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = TRUE)) %>%
  ungroup() %>%
  arrange(study_arm, beta_quartile)

write.table(beta_table_v2, file = paste0(datestring, "_MDCF_POC_participant_beta_WLZ_quartiles.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

rm(beta_table)

beta_table_v2 %>%
  group_by(study_arm, beta_quartile) %>%
  summarize(mean_beta_wlz = mean(coef))

MDCF2_Q4 <- beta_table_v2 %>%
  filter(study_arm == "X" & beta_quartile == "Q4")

MDCF2_Q3 <- beta_table_v2 %>%
  filter(study_arm == "X" & beta_quartile == "Q3")

ggplot(map_merge %>% filter(PID %in% MDCF2_Q3$PID), aes(x = study_week, y = wlz, group = PID)) +
  geom_line() +
  facet_wrap(~ PID)
