################################################################################
#
# Contamination removal benchmarking - SourceTracker 
#
################################################################################

BM_sourcetracker <- function(file_directory, otus, metadata, Mock_info,
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
  
  # Tune alpha parameters - not working with only one given source environment
  #tune_test_res <- my_tune.st(otus = otus, envs = envs,
  #                            individual.samples = FALSE)

  # Define SourceTracker training data
  set.seed(1)
  st_train <- sourcetracker(otus[train.ix, ], envs[train.ix], 
                            rarefaction_depth = raref_depth)
  
  # Increase alpha parameters to avoid overfitting
  # Scenario 1 - with default parameters
  set.seed(1)
  st_res_default <- predict.sourcetracker(
    stobj = st_train, test = otus[test.ix, ], rarefaction_depth = raref_depth,
    alpha1 = 0.001, alpha2 = 0.1, beta = 10, full.results = TRUE)
  
  # Scenario 2 - alpha1 = 0.1
  set.seed(1)
  st_res_alpha1_0.1 <- predict.sourcetracker(
    stobj = st_train, test = otus[test.ix, ], rarefaction_depth = raref_depth,
    alpha1 = 0.1, alpha2 = 0.1, beta = 10, full.results = TRUE)

  # Scenario 3 - alpha2 = 0.001
  set.seed(1)
  st_res_alpha2_0.001 <- predict.sourcetracker(
    stobj = st_train, test = otus[test.ix, ], rarefaction_depth = raref_depth,
    alpha1 = 0.001, alpha2 = 0.001, beta = 10, full.results = TRUE)

  # Scenario 4 - beta = 0.001
  set.seed(1)
  st_res_beta_0.001 <- predict.sourcetracker(
    stobj = st_train, test = otus[test.ix, ], rarefaction_depth = raref_depth,
    alpha1 = 0.001, alpha2 = 0.1, beta = 0.001, full.results = TRUE)

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
    list(st_res_default, st_res_alpha1_0.1, st_res_alpha2_0.001, 
         st_res_beta_0.001),
    nm = c("default", "alpha1 = 0.1", "alpha2 = 0.001", "beta = 0.001"))
  
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
