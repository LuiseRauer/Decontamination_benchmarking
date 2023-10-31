# Decontamination benchmarking

The code provided here benchmarks six bioinformatic decontamination algorithms for microbiome sequencing data in three mock datasets. This work is published as "Benchmarking MicrobIEM – a user-friendly tool for decontamination of microbiome sequencing data" in BMC Biology (details follow).  

The tested algorithms include the established frequency filter, Decontam frequency filter, Decontam prevalence filter, SourceTracker, presence filter, and the novel ratio and span filter implemented in our tool MicrobIEM (https://github.com/LuiseRauer/MicrobIEM).  

R code files can be executed in any order. Files starting with "BM_" contain the functions for running the respective decontamination algorithms and evaluating the number of true positives, true negatives, false positives and false negatives for any given data set and gold standard classification. Benchmarking of the decontamination algorithms for each dataset can be found in scripts starting with "1-3".  

- 0_Staggered-mock-A_DADA2.R: DADA2-processing of raw FASTQ files of the staggered mock A  
- 1_Even-mock_benchmarking.R: application of decontamination algorithms to the even mock community, as provided by Karstens et al.  
- 2_Staggered-mock-A_benchmarking.R: application of decontamination algorithms to the staggered mock community A, as pre-processed in script 0  
- 3_Staggered-mock-B_benchmarking.R: application of decontamination algorithms to the staggered mock community B, as provided by Rauer & De Tomassi et al.  
