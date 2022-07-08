################################################################################
#
# Contamination removal benchmarking - Decontam prevalence filter 
#
################################################################################

BM_decontam_prev <- function(otus_rel, metadata, Mock_info, per_dilution, dmax) {

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
  # Evaluate results for contamination removal - Prevalence filter
  # ------------------------------------------------------------------------------
  
  # Define thresholds
  thresholds_prev <- c(0.2, 0.4, 0.6, 0.8, 1)
  
  # Create an empty dataframe for results
  res_DecontamPrev <- 
    setNames(data.frame(matrix(ncol = 6, nrow = 0)), 
             c("Sample", "TP", "TN", "FP", "FN", "Filter"))
  
  # loop for prevalence method, reseq dataset
  for(i in thresholds_prev) {
    #i <- thresholds_prev[1]
    for(k in 0:dmax){ 
      #k <- 0
      k_samples <- metadata[grepl(paste0("D", k), metadata$Dilution), "Sample_ID"]
      # Decide if analysis should be performed on each dilution individually 
      if (per_dilution) {
        otus_rel_used <- as.matrix(otus_rel[rownames(otus_rel) %in% c(k_samples, Control_vec), ])
        metadata_used <- metadata[metadata$Sample_ID %in% c(k_samples, Control_vec), ]
      } else {
        otus_rel_used <- as.matrix(otus_rel)
        metadata_used <- metadata 
      }  
      is_control <- ifelse(metadata_used$Sample_type == "SAMPLE", FALSE, TRUE)
      contamdf <- isContaminant(seqtab = otus_rel_used,  
                                method = "prevalence", neg = is_control, 
                                threshold = i)
      contamdf <- merge(t(otus_rel_used)[, k_samples, drop = FALSE], # !!!
                        contamdf[, "contaminant", drop = FALSE], 
                        by = 0, all.x = TRUE)
      contamdf <- merge(contamdf, Mock_info,
                        by.x = "Row.names", by.y = 0, all.x = TRUE)
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
      if(identical(TN, numeric(0))) {TN = rep(0, length(k_samples))} # !!!
      if(identical(FN, numeric(0))) {FN = rep(0, length(k_samples))} # !!!
      if(identical(TP, numeric(0))) {TP = rep(0, length(k_samples))} # !!!
      if(identical(FP, numeric(0))) {FP = rep(0, length(k_samples))} # !!!
      res_DecontamPrev <- 
        rbind.data.frame(res_DecontamPrev,
                         data.frame(Sample = k_samples,
                                    TP = TP, TN = TN, FP = FP, FN = FN,
                                    Filter = paste0("Decontam (prevalence), thr = ", i)))
    }
  }
  return(res_DecontamPrev)
}
