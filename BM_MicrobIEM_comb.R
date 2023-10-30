################################################################################
#
# Contamination removal benchmarking - MicrobIEM span filter 
#
################################################################################

BM_microbiem_comb <- function(otus_rel, metadata, Mock_info, dmax) {
  
  # ------------------------------------------------------------------------------
  # Set up the workspace
  # ------------------------------------------------------------------------------
  
  # Load required packages
  library(tidyverse)
  # Download MicrobIEM: https://github.com/LuiseRauer/MicrobIEM/releases/tag/v0.7
  # Unzip the folder and load the script
  source(paste0(file_directory, 
                "Input/MicrobIEM/MicrobIEM-0.7/MicrobIEM_decontamination.r"))
  
  # ----------------------------------------------------------------------------
  # Prepare the data
  # ----------------------------------------------------------------------------
  
  # Define a vector of negative controls
  Control_vec <- metadata[metadata$Sample_type == "NEG2", "Sample_ID"]
  
  # ------------------------------------------------------------------------------
  # Evaluate results for contamination removal
  # ------------------------------------------------------------------------------
  
  # Create an empty dataframe for results
  res_MicrobiemComb <- 
    setNames(data.frame(matrix(ncol = 6, nrow = 0)), 
             c("Sample", "TP", "TN", "FP", "FN", "Filter"))
  
  ratio_thresholds <- c(2, 1.5, 1, 0.5, 0.1, NA)
  span_thresholds <- c(NA, seq(0, 1, 1/length(Control_vec))[2:(length(Control_vec)+1)])
  
  for (i in ratio_thresholds) {
    for (j in span_thresholds) {
      #i = 5
      # Select only one dilution
      for (k in 0:dmax) { 
        #k = 4
        # Subset dilution
        k_samples <- metadata[grepl(paste0("D", k), metadata$Dilution), "Sample_ID"]
        contamdf <- MicrobIEM_decontamination(
          feature_table = as.data.frame(t(otus_rel)), SAMPLE = k_samples, 
          NEG2 = Control_vec, ratio_NEG2_threshold = i, span_NEG2_threshold = j)
        contamdf <- merge(t(otus_rel)[, k_samples, drop = FALSE], # !!!
                          contamdf[, "is_contaminant", drop = FALSE], 
                          by = 0, all.x = TRUE)
        contamdf <- merge(contamdf, Mock_info,
                          by.x = "Row.names", by.y = 0, all.x = TRUE)
        # Calculation of accuracy, sensitivity, specificity
        # Sum OTUs classified as mock/cont. per environment to contingency table
        otus_rel_sum <- contamdf %>% 
          group_by(is_contaminant, Contaminant) %>% 
          summarise_at(.vars = all_of(k_samples), .funs = sum) %>% as.data.frame()
        # Define TN/FN/TP/FP
        TN = filter(otus_rel_sum, is_contaminant == "FALSE" & Contaminant == "Mock") %>%
          select(all_of(k_samples)) %>% unlist()
        FN = filter(otus_rel_sum, is_contaminant == "FALSE" & Contaminant == "Contaminant") %>% 
          select(all_of(k_samples)) %>% unlist()
        TP = filter(otus_rel_sum, is_contaminant == "TRUE" & Contaminant == "Contaminant") %>%
          select(all_of(k_samples)) %>% unlist()
        FP = filter(otus_rel_sum, is_contaminant == "TRUE" & Contaminant == "Mock") %>%
          select(all_of(k_samples)) %>% unlist()
        # Add missing categories
        if(all(is.na(TN))) {TN = rep(0, length(k_samples))} 
        if(all(is.na(FN))) {FN = rep(0, length(k_samples))} 
        if(all(is.na(TP))) {TP = rep(0, length(k_samples))} 
        if(all(is.na(FP))) {FP = rep(0, length(k_samples))} 
        res_MicrobiemComb <- 
          rbind.data.frame(res_MicrobiemComb,
                           data.frame(Sample = k_samples,
                                      TP = TP, TN = TN, FP = FP, FN = FN,
                                      Filter = paste0("MicrobIEM, ratio = ", i,
                                                      "; span = ", j)))
      }
    }
  }
  return(res_MicrobiemComb)
}
