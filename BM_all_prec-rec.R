################################################################################
#
# Contamination removal benchmarking - SourceTracker 
#
################################################################################

BM_sourcetracker_PR <- function(file_directory, otus, metadata, Mock_info,
                                raref_depth) {
  
  # ------------------------------------------------------------------------------
  # Set up the workspace
  # ------------------------------------------------------------------------------
  
  # Load required packages
  library(tidyverse)
  # Download SourceTracker: https://github.com/danknights/sourcetracker/releases/tag/v1.0.1
  # Unzip the folder and load the script
  source(paste0(file_directory, 
                "Input/SourceTracker/sourcetracker-1.0.1/src/SourceTracker.r"))
  
  # ------------------------------------------------------------------------------
  # Define source samples as training data, sink samples as test data
  # ------------------------------------------------------------------------------
  
  # Extract the source environments and source/sink indices
  train.ix <- which(metadata$SourceSink == "Source")
  test.ix <- which(metadata$SourceSink == "Sink")
  envs <- as.factor(metadata$Sample_type)
  
  # ------------------------------------------------------------------------------
  # Run SourceTracker algorithm - estimate source proportions in test data
  # ------------------------------------------------------------------------------
  print("E1")
  # Tune alpha parameters - not working with only one given source environment
  #tune_test_res <- my_tune.st(otus = otus, envs = envs,
  #                            individual.samples = FALSE)
  
  # Define SourceTracker training data
  set.seed(1)
  st_train <- sourcetracker(otus[train.ix, ], envs[train.ix], 
                            rarefaction_depth = raref_depth)
  
  thresholds <- format(c(10^seq(-3, 2, 1)), scientific = T)
  
  # Increase alpha parameters to avoid overfitting
  # Scenario 1 - with default parameters
  for (i in thresholds) {
    set.seed(1)
    st_res_a1 <- predict.sourcetracker(
      stobj = st_train, test = otus[test.ix, ], rarefaction_depth = raref_depth,
      alpha1 = as.numeric(i), alpha2 = 0.1, beta = 10, full.results = TRUE)
    assign(paste0("st_res_a1_", i), st_res_a1)
  }
  
  # Scenario 2 - alpha1 = 0.1
  for (i in thresholds) {
    set.seed(1)
    st_res_a2 <- predict.sourcetracker(
      stobj = st_train, test = otus[test.ix, ], rarefaction_depth = raref_depth,
      alpha1 = 0.1, alpha2 = as.numeric(i), beta = 10, full.results = TRUE)
    assign(paste0("st_res_a2_", i), st_res_a2)
  }  
  
  # Scenario 3 - alpha2 = 0.001
  for (i in thresholds) {
    set.seed(1)
    st_res_b1 <- predict.sourcetracker(
      stobj = st_train, test = otus[test.ix, ], rarefaction_depth = raref_depth,
      alpha1 = 0.001, alpha2 = 0.1, beta = as.numeric(i), full.results = TRUE)
    assign(paste0("st_res_b1_", i), st_res_b1)
  }  
  # Scenario 4 - beta = 0.001
  #set.seed(1)
  #st_res_beta_0.001 <- predict.sourcetracker(
  #  stobj = st_train, test = otus[test.ix, ], rarefaction_depth = raref_depth,
  #  alpha1 = 0.001, alpha2 = 0.1, beta = 0.001, full.results = TRUE)
  
  # Structure of full.results:
  ###dim(st_res_default$full.results)
  # 10 = total number of draws
  # 2 = number of source envs
  # 1384 = number of taxa/OTUs/ASVs
  # 9 = number of sink samples
  
  # ------------------------------------------------------------------------------
  # Evaluate results for contamination removal
  # ------------------------------------------------------------------------------
  
  # Create a list of SourceTracker results
  st_list <- set_names(
    mget(c(paste0("st_res_a1_", thresholds), paste0("st_res_a2_", thresholds),
           paste0("st_res_b1_", thresholds))),
    nm = paste0(rep(c("a1 = ", "a2 = ", "b1 = "), each = length(thresholds)), 
                rep(thresholds, 3)))
  
  # Define empty dataframe to add results in following loop
  res_SourceTracker <- 
    setNames(data.frame(matrix(ncol = 6, nrow = 0)), 
             c("Sample", "TP", "TN", "FP", "FN", "Filter"))
  
  for (k in 1:length(st_list)) {
    #k = 1
    j <- st_list[[k]]
    # Add proper dimension names to SourceTracker output
    st_res <- array(
      j$full.results, 
      dim = dim(j$full.results), 
      dimnames = list(
        paste0("draw_", 1:dim(j$full.results)[1]), # Draws
        j$train.envs, # Environments
        colnames(otus), # Taxa/OTUs/ASVs
        j$samplenames)) # Samples
    # Get average counts from all draws of SourceTracker output
    st_res_mean <- apply(st_res, c(2, 3, 4), mean) 
    # Summarise results per sample
    for (i in 1:dim(st_res_mean)[3]) { #18
      #i = 1
      # Merge SourceTracker result per sample with gold standard classification
      st_res_mean_class <- merge(
        t(st_res_mean[, , i]), Mock_info, by = 0, all.x = TRUE)
      # Calculation of accuracy, sensitivity, specificity
      # Sum OTUs classified as mock/cont. per environment to contingency table
      st_res_mean_sum <- st_res_mean_class %>% 
        group_by(Contaminant) %>% 
        summarise_at(.vars = c("NEG2", "Unknown"), .funs = sum) %>% as.data.frame()
      st_res_mean_sum <- column_to_rownames(st_res_mean_sum, "Contaminant")
      TP = st_res_mean_sum["Contaminant", "NEG2"]
      TN = st_res_mean_sum["Mock", "Unknown"]
      FP = st_res_mean_sum["Mock", "NEG2"]
      FN = st_res_mean_sum["Contaminant", "Unknown"]
      res_SourceTracker <- rbind.data.frame(
        res_SourceTracker,
        data.frame(Sample = dimnames(st_res_mean)[[3]][i], # Sample name
                   TP = TP, TN = TN, FP = FP, FN = FN, 
                   Filter = paste0("SourceTracker, ", names(st_list)[k]))
      ) 
    }
  }
  return(res_SourceTracker)
}

################################################################################
#
# Contamination removal benchmarking - Decontam frequency filter 
#
################################################################################
library(tidyverse)
BM_decontam_freq_PR <- function(otus_rel, metadata, Mock_info, per_dilution, dmax) {
  
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
  thresholds_freq <- seq(0, 1, 0.05)
  
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

################################################################################
#
# Contamination removal benchmarking - Decontam prevalence filter 
#
################################################################################

BM_decontam_prev_PR <- function(otus_rel, metadata, Mock_info, per_dilution, dmax) {
  
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
  thresholds_prev <- seq(0, 1, 0.05)
  
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

################################################################################
#
# Contamination removal benchmarking - MicrobIEM ratio filter 
#
################################################################################

BM_microbiem_ratio_PR <- function(otus_rel, metadata, Mock_info, dmax) {
  
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
  res_MicrobiemRatio <- 
    setNames(data.frame(matrix(ncol = 6, nrow = 0)), 
             c("Sample", "TP", "TN", "FP", "FN", "Filter"))
  
  thresholds <- sort(c(10^seq(-3, 3, 1), 2*10^seq(-3, 3, 1), 5*10^seq(-3, 3, 1)))
  
  for (i in thresholds) {
    #i = 5
    # Select only one dilution
    for (k in 0:dmax) { 
      print(paste(k, i, sep = "___"))
      #k = 4
      # Subset dilution
      k_samples <- metadata[grepl(paste0("D", k), metadata$Dilution), "Sample_ID"]
      contamdf <- MicrobIEM_decontamination(
        feature_table = as.data.frame(t(otus_rel)), SAMPLE = k_samples, 
        NEG2 = Control_vec, ratio_NEG2_threshold = i, span_NEG2_threshold = NA)
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
      res_MicrobiemRatio <- 
        rbind.data.frame(res_MicrobiemRatio,
                         data.frame(Sample = k_samples,
                                    TP = TP, TN = TN, FP = FP, FN = FN,
                                    Filter = paste0("MicrobIEM, ratio = ", i)))
    }
  }
  return(res_MicrobiemRatio)
}

################################################################################
#
# Contamination removal benchmarking - MicrobIEM span filter 
#
################################################################################

BM_microbiem_span_PR <- function(otus_rel, metadata, Mock_info, dmax) {
  
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
  res_MicrobiemSpan <- 
    setNames(data.frame(matrix(ncol = 6, nrow = 0)), 
             c("Sample", "TP", "TN", "FP", "FN", "Filter"))
  
  thresholds <- seq(0, 1, 1/length(Control_vec))[2:(length(Control_vec)+1)]
  
  for (i in thresholds) {
    #i = 5
    # Select only one dilution
    for (k in 0:dmax) { 
      #k = 4
      # Subset dilution
      k_samples <- metadata[grepl(paste0("D", k), metadata$Dilution), "Sample_ID"]
      contamdf <- MicrobIEM_decontamination(
        feature_table = as.data.frame(t(otus_rel)), SAMPLE = k_samples, 
        NEG2 = Control_vec, ratio_NEG2_threshold = NA, span_NEG2_threshold = i)
      contamdf <- merge(t(otus_rel)[, k_samples, drop = FALSE], 
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
      res_MicrobiemSpan <- 
        rbind.data.frame(res_MicrobiemSpan,
                         data.frame(Sample = k_samples,
                                    TP = TP, TN = TN, FP = FP, FN = FN,
                                    Filter = paste0("MicrobIEM, span = ", i)))
    }
  }
  return(res_MicrobiemSpan)
}

################################################################################
#
# Contamination removal benchmarking - Frequency filter 
#
################################################################################

BM_frequency_PR <- function(otus_rel, metadata, Mock_info, dmax) {
  
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
  thresholds <- sort(c(10^seq(-1, -8, -1), 5*10^seq(-1, -8, -1)))
  
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
