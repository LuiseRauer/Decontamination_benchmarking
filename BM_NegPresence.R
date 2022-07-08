################################################################################
#
# Contamination removal benchmarking - Negative control filter 
#
################################################################################

BM_negpresence <- function(otus_rel, metadata, Mock_info, dmax) {
  
  # ------------------------------------------------------------------------------
  # Set up the workspace
  # ------------------------------------------------------------------------------
  
  # Load required packages
  library(tidyverse)
  
  # ------------------------------------------------------------------------------
  # Prepare the data
  # ------------------------------------------------------------------------------

  # Define a vector of negative controls
  Control_vec <- metadata[metadata$Sample_type == "NEG2", "Sample_ID"]
  # Add mock/contaminant status by NEG2 presence/absence
  otus_NEG2 <- as.data.frame(t(otus_rel))
  otus_NEG2["Filter"] <- ifelse(
    rowSums(otus_NEG2[, Control_vec, drop = FALSE]) > 0, "Cont_f", "Mock_f")# !!!
  # Add mock/contaminant status by gold standard
  otus_NEG2 <- merge(otus_NEG2, Mock_info, by = 0, all.x = TRUE)
  otus_NEG2 <- column_to_rownames(otus_NEG2, "Row.names")
  
  #otus_NEG2["Contaminant"] <- ifelse(row.names(otus_NEG2) %in% mock_taxa,
  #                                   "Mock", "Contaminant")
  
  # ------------------------------------------------------------------------------
  # Evaluate results for contamination removal
  # ------------------------------------------------------------------------------
  
  # Define empty vectors to add results in following loop
  res_NegPresence <- 
    setNames(data.frame(matrix(ncol = 6, nrow = 0)), 
             c("Sample", "TP", "TN", "FP", "FN", "Filter"))
  for (k in 0:dmax) {
    #k = 0
    # Subset dilution
    k_samples <- metadata[grepl(paste0("D", k), metadata$Dilution), "Sample_ID"]
    # Calculation of accuracy, sensitivity, specificity
    # Sum reads classified as mock/cont. per environment to contingency table
    otus_NEG2_sum <- otus_NEG2 %>% 
      group_by(Filter, Contaminant) %>% 
      summarise_at(.vars = all_of(k_samples), .funs = sum) %>% as.data.frame()
    # Define TN/FN/TP/FP
    TN = filter(otus_NEG2_sum, Filter == "Mock_f" & Contaminant == "Mock") %>%
      select(all_of(k_samples)) %>% unlist()
    FN = filter(otus_NEG2_sum, Filter == "Mock_f" & Contaminant == "Contaminant") %>%
      select(all_of(k_samples)) %>% unlist()
    TP = filter(otus_NEG2_sum, Filter == "Cont_f" & Contaminant == "Contaminant") %>%
      select(all_of(k_samples)) %>% unlist()
    FP = filter(otus_NEG2_sum, Filter == "Cont_f" & Contaminant == "Mock") %>%
      select(all_of(k_samples)) %>% unlist
    # Add missing categories
    if(all(is.na(TN))) {TN = rep(0, length(k_samples))} 
    if(all(is.na(FN))) {FN = rep(0, length(k_samples))} 
    if(all(is.na(TP))) {TP = rep(0, length(k_samples))} 
    if(all(is.na(FP))) {FP = rep(0, length(k_samples))} 
    res_NegPresence <- 
      rbind.data.frame(res_NegPresence,
                       data.frame(Sample = k_samples,
                                  TP = TP, TN = TN, FP = FP, FN = FN, 
                                  Filter = "Presence in NEG2",
                                  row.names = NULL))
  }
  return(res_NegPresence)
}
