################################################################################
#
# Contamination removal benchmarking - Decontam frequency filter 
#
################################################################################

BM_decontam_freq <- function(otus_rel, metadata, Mock_info, per_dilution, dmax) {
  
  # ------------------------------------------------------------------------------
  # Set up the workspace
  # ------------------------------------------------------------------------------
  
  # Installation (https://bioconductor.org/packages/release/bioc/html/decontam.html)
  ###if (!requireNamespace("BiocManager", quietly = TRUE))
  ###  install.packages("BiocManager")
  ###BiocManager::install("decontam")
  
  # Load required packages
  library(tidyverse)
  library(decontam); packageVersion("decontam")
  
  # ------------------------------------------------------------------------------
  # Prepare the data
  # ------------------------------------------------------------------------------
  
  # Define a vector of negative controls
  Control_vec <- metadata[metadata$Sample_type == "NEG2", "Sample_ID"]
  
  # ------------------------------------------------------------------------------
  # Evaluate results for contamination removal - Frequency filter
  # ------------------------------------------------------------------------------
  
  # Define thresholds
  thresholds_freq <- c(0.2, 0.4, 0.6, 0.8, 1)
  
  # Create an empty dataframe for results
  res_DecontamFreq <- 
    setNames(data.frame(matrix(ncol = 6, nrow = 0)), 
             c("Sample", "TP", "TN", "FP", "FN", "Filter"))
  
  for(i in thresholds_freq) {
    #i <- thresholds_freq[1]
    for(k in 0:dmax){ 
      #k <- 6
      k_samples <- metadata[grepl(paste0("D", k), metadata$Dilution), "Sample_ID"]
      # Decide if analysis should be performed on each dilution individually 
      if (per_dilution) {
        otus_rel_used <- as.matrix(otus_rel[rownames(otus_rel) %in% k_samples, ])
        metadata_used <- metadata[metadata$Sample_ID %in% k_samples, ]
      } else {
        otus_rel_used <- as.matrix(otus_rel[!rownames(otus_rel) %in% Control_vec, ])
        metadata_used <- metadata[!metadata$Sample_ID %in% Control_vec, ]
      }
      quant_reading <- metadata_used$DNA_conc 
      #is_control <- ifelse(metadata_used$Sample_type == "SAMPLE", FALSE, TRUE)
      contamdf <- isContaminant(seqtab = otus_rel_used, conc = quant_reading, 
                                method = "frequency", #neg = is_control, 
                                threshold = i)
      contamdf <- merge(t(otus_rel_used)[, k_samples, drop = FALSE], # !!!
                        contamdf[, "contaminant", drop = FALSE], 
                        by = 0, all.x = TRUE)
      contamdf <- merge(contamdf, Mock_info, by.x = "Row.names", by.y = 0, 
                        all.x = TRUE)
      # With sum for cont/mock reduction
      contamdf_sum <- contamdf %>% 
        group_by(contaminant, Contaminant) %>% 
        summarise_at(.vars = all_of(k_samples), .funs = sum) %>% as.data.frame()
      TP <- filter(contamdf_sum, contaminant == "TRUE", Contaminant == "Contaminant") %>%
        select(all_of(k_samples)) %>% unlist
      TN <- filter(contamdf_sum, contaminant == "FALSE", Contaminant == "Mock") %>%
        select(all_of(k_samples)) %>% unlist
      FP <- filter(contamdf_sum, contaminant == "TRUE", Contaminant == "Mock") %>%
        select(all_of(k_samples)) %>% unlist
      FN <- filter(contamdf_sum, contaminant == "FALSE", Contaminant == "Contaminant") %>%
        select(all_of(k_samples)) %>% unlist
      # Add missing categories
      if(identical(TN, numeric(0))) {TN = rep(0, length(k_samples))} 
      if(identical(FN, numeric(0))) {FN = rep(0, length(k_samples))} 
      if(identical(TP, numeric(0))) {TP = rep(0, length(k_samples))} 
      if(identical(FP, numeric(0))) {FP = rep(0, length(k_samples))} 
      res_DecontamFreq <- 
        rbind.data.frame(res_DecontamFreq,
                         data.frame(Sample = k_samples,
                                    TP = TP, TN = TN, FP = FP, FN = FN,
                                    Filter = paste0("Decontam (frequency), thr = ", i)))
    }
  }
  return(res_DecontamFreq)
}
