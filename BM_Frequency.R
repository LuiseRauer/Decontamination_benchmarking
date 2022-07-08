################################################################################
#
# Contamination removal benchmarking - Frequency filter 
#
################################################################################

BM_frequency <- function(otus_rel, metadata, Mock_info, dmax) {
  
  # ------------------------------------------------------------------------------
  # Set up the workspace
  # ------------------------------------------------------------------------------
  
  # Load required packages
  library(tidyverse)
  
  # ------------------------------------------------------------------------------
  # Evaluate results for contamination removal
  # ------------------------------------------------------------------------------
  
  # Create an empty dataframe for results
  res_FrequencyFilter <- 
    setNames(data.frame(matrix(ncol = 6, nrow = 0)), 
             c("Sample", "TP", "TN", "FP", "FN", "Filter"))
  
  # Define thresholds for frequency filter
  thresholds <- c(0.0001, 0.001, 0.01)
  
  for (i in thresholds) {
    #i = 0.0001
    #print(i)
    # Summarise results per dilution
    for (k in 0:dmax) { # selection of k automatically removes NEG2 samples # !!!
      #k = 4
      #print(k)
      # Subset dilution
      k_samples <- metadata[grepl(paste0("D", k), metadata$Dilution), "Sample_ID"]
      otus_rel_class <- as.data.frame(t(otus_rel[k_samples, ]))
      # Apply frequency filter
      otus_rel_class$Filter <- ifelse(apply(otus_rel_class, 1, max) >= i,
                                      "Mock_f", "Cont_f")
      # Merge filter result per filter i with mock/contaminant classification
      otus_rel_class <- merge(otus_rel_class, Mock_info, by = 0, all.x = TRUE)
      # Calculation of accuracy, sensitivity, specificity
      # Sum OTUs classified as mock/cont. per environment to contingency table
      otus_rel_sum <- otus_rel_class %>% 
        group_by(Filter, Contaminant) %>% 
        summarise_at(.vars = all_of(k_samples), .funs = sum) %>% as.data.frame()
      # Define TN/FN/TP/FP
      TN = filter(otus_rel_sum, Filter == "Mock_f" & Contaminant == "Mock") %>%
        select(all_of(k_samples)) %>% unlist()
      FN = filter(otus_rel_sum, Filter == "Mock_f" & Contaminant == "Contaminant") %>% 
        select(all_of(k_samples)) %>% unlist()
      TP = filter(otus_rel_sum, Filter == "Cont_f" & Contaminant == "Contaminant") %>%
        select(all_of(k_samples)) %>% unlist()
      FP = filter(otus_rel_sum, Filter == "Cont_f" & Contaminant == "Mock") %>%
        select(all_of(k_samples)) %>% unlist()
      # Add missing categories
      if(all(is.na(TN))) {TN = rep(0, length(k_samples))} 
      if(all(is.na(FN))) {FN = rep(0, length(k_samples))} 
      if(all(is.na(TP))) {TP = rep(0, length(k_samples))} 
      if(all(is.na(FP))) {FP = rep(0, length(k_samples))} 
      res_FrequencyFilter <- 
        rbind.data.frame(res_FrequencyFilter,
                         data.frame(Sample = k_samples,
                                    TP = TP, TN = TN, FP = FP, FN = FN, 
                                    Filter = paste0("Frequency filter, thr = ", 
                                                    format(i, scientific = FALSE)),
                                    row.names = NULL))
    }
  }
  return(res_FrequencyFilter)
}
