# functions

# Some samples are sequenced across multiple runs, so aggregate their data here
#map.merge <- map[!duplicated(map$new_SampleID),]
# map.merge <- map %>%
#   group_by(new_SampleID) %>%
#   summarise(X.SampleID = ifelse(length(unique(X.SampleID)) > 1, paste(X.SampleID, collapse = ";"), unique(as.character(X.SampleID))),
#             BarcodeSequence = ifelse(length(unique(BarcodeSequence)) > 1, paste(BarcodeSequence, collapse = ";"), unique(as.character(BarcodeSequence))),
#             LinkerPrimerSequence = ifelse(length(unique(LinkerPrimerSequence)) > 1, paste(LinkerPrimerSequence, collapse = ";"), unique(as.character(LinkerPrimerSequence))),
#             #            new_SampleID = paste(new_SampleID, collapse = ";"),
#             old_SampleID = ifelse(length(unique(old_SampleID)) > 1, paste(old_SampleID, collapse = ";"), unique(as.character(old_SampleID))),
#             cohort_1 = ifelse(length(unique(cohort_1)) > 1, paste(cohort_1, collapse = ";"), unique(as.character(cohort_1))),
#             participant = ifelse(length(unique(participant)) > 1, paste(participant, collapse = ";"), unique(as.character(participant))),
#             donor_age = ifelse(length(unique(donor_age)) > 1, paste(donor_age, collapse = ";"), unique(as.character(donor_age))),
#             sex = ifelse(length(unique(sex)) > 1, paste(sex, collapse = ";"), unique(as.character(sex))),
#             whz = ifelse(length(unique(whz)) > 1, paste(whz, collapse = ";"), unique(as.character(whz))),
#             waz = ifelse(length(unique(waz)) > 1, paste(waz, collapse = ";"), unique(as.character(waz))),
#             haz = ifelse(length(unique(haz)) > 1, paste(haz, collapse = ";"), unique(as.character(haz))),
#             health_status = ifelse(length(unique(health_status)) > 1, paste(health_status, collapse = ";"), unique(as.character(health_status))),
#             expected_reads = ifelse(length(unique(expected_reads)) > 1, paste(expected_reads, collapse = ";"), unique(as.character(expected_reads))),
#             st = ifelse(length(unique(st)) > 1, paste(st, collapse = ";"), unique(as.character(st))),
#             dada2 = ifelse(length(unique(dada2)) > 1, paste(dada2, collapse = ";"), unique(as.character(dada2))),
#             seq_scheme = ifelse(length(unique(seq_scheme)) > 1, paste(seq_scheme, collapse = ";"), unique(as.character(seq_scheme))),
#             Description = ifelse(length(unique(Description)) > 1, paste(Description, collapse = ";"), unique(as.character(Description)))) %>%
#   as.data.frame(stringsAsFactors = F) %>%
#   sample_data()
# row.names(map.merge) <- map.merge$new_SampleID
# ps1 <- merge_samples(ps0, "new_SampleID", fun = sum)
# sample_data(ps1) <- sample_data(map.merge)
# nsamples(ps1)


get_coef_p <- function(model, term) {
  coef_df <- as.data.frame(summary(model)$coefficients)
  row.names(coef_df) <- unlist(dimnames(coef_df)[1])
  colnames(coef_df) <- unlist(dimnames(coef_df)[2])
  coef <- coef_df[term, "Estimate"]
  pval <- coef_df[term, "Pr(>|t|)"]
  return(data.frame(term = term, coef = coef, pval = pval))
}

#lm.beta <- 
#  function (MOD) 
#  {
#    b <- summary(MOD)$coef[-1, 1]
#    sx <- sd(MOD$model[-1])
#    sy <- sd(MOD$model[1])
#    beta <- b * sx/sy
#    return(beta)
#  }

get_anova_p <- function(anova, term) {
  anova_df <- as.data.frame(anova) %>%
    rownames_to_column(var = "term")
  pval <- anova_df[anova_df$term == term, "Pr(>F)"]
  return(pval)
}

get_confint <- function(model, term) {
  confint_df <- as.data.frame(confint(profile(model), level = 0.95)) %>%
    rownames_to_column("term")
  lci <- confint_df[confint_df$term == term,]$`2.5 %`
  uci <- confint_df[confint_df$term == term,]$`97.5 %`
  if (length(lci) > 0) {
    return(paste(lci, uci, sep = ";"))
  } else {
    return(paste(NA, NA, sep = ";"))
  }
  
}

get_coef <- function(model, term) {
  coef_df <- as.data.frame(summary(model)$coefficients)
  row.names(coef_df) <- unlist(dimnames(coef_df)[1])
  colnames(coef_df) <- unlist(dimnames(coef_df)[2])
  coef <- coef_df[term, "Estimate"]
  #  pval <- coef_df[term, "Pr(>|t|)"]
  if (is.na(coef)) {
    print(paste0("WARNING: term ", term, " does not exist in model!"))
    return(NA)
  } else {
    return(coef)
  }
}

get_std_coef <- function(model, term, response, predictor) {
  coef_df <- as.data.frame(summary(model)$coefficients)
  row.names(coef_df) <- unlist(dimnames(coef_df)[1])
  colnames(coef_df) <- unlist(dimnames(coef_df)[2])
  coef <- coef_df[term, "Estimate"]
  coef_std <- coef*(sd(predictor)/sd(response))
  #  pval <- coef_df[term, "Pr(>|t|)"]
  if (is.na(coef_std)) {
    print(paste0("WARNING: term ", term, " does not exist in model!"))
    return(NA)
  } else {
    return(coef_std)
  }
}

get_intercept <- function(model) {
  coef_df <- as.data.frame(summary(model)$coefficients)
  row.names(coef_df) <- unlist(dimnames(coef_df)[1])
  colnames(coef_df) <- unlist(dimnames(coef_df)[2])
  coef <- coef_df["(Intercept)", "Estimate"]
  #  pval <- coef_df[term, "Pr(>|t|)"]
  if (is.na(coef)) {
    stop(paste0("ERROR: term ", term, " does not exist in model!"))
  } else {
    return(coef)
  }
}

sem <- function(x, na.rm) sd(x, na.rm = na.rm)/sqrt(length(x[!is.na(x)]))

get_model_se <- function(model, term) {
  coef_df <- as.data.frame(summary(model)$coefficients)
  row.names(coef_df) <- unlist(dimnames(coef_df)[1])
  colnames(coef_df) <- unlist(dimnames(coef_df)[2])
  coef <- coef_df[term, "Std. Error"]
  if (is.na(coef)) {
    print(paste0("WARNING: term ", term, " does not exist in model!"))
    return(NA)
  } else {
    return(coef)
  }
}

#adapted from: https://stackoverflow.com/questions/17171148/non-redundant-version-of-expand-grid
expand.grid.unique <- function(x, y, include.equals = FALSE) {
  x <- unique(x)
  y <- unique(y)
  g <- function(i) {
    z <- setdiff(y, x[seq_len(i-include.equals)] )
    if (length(z) ) {
      cbind(x[i], z, deparse.level=0)
    }
  }
  do.call(rbind, lapply(seq_along(x), g))
}

log10_plus_pseudo <- function(x) {log10(x + 0.5)}
