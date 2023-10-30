################################################################################
#
# Contamination removal benchmarking - Identifying correct ASVs with Sanger data 
#
################################################################################

# Define directory
file_directory <- "C:/Users/rauerlui/PhD/Projects/2020-03_MicrobIEM/2022.07 Decontamination benchmarking/Decontamination_benchmarking/"

# Read ASV table
Dmock_ASVtable <- read.csv(
  paste0(file_directory, "Input/DilutionMock/Dmock_ASVtable.csv"), sep = ";",
         check.names = FALSE)
# Sort ASVs by ASV_ID
Dmock_ASVtable <- Dmock_ASVtable[order(Dmock_ASVtable$OTU_ID), ]
row.names(Dmock_ASVtable) <- NULL

# ------------------------------------------------------------------------------
# Set up the workspace
# ------------------------------------------------------------------------------

# Installation
###if (!require("BiocManager", quietly = TRUE))
###  install.packages("BiocManager")
###BiocManager::install("sangerseqR")

library(sangerseqR)
library(stringr)

# Path of input FASTQ files (must be demultiplexed, can be fastq.gz)
path_sanger <- paste0(file_directory, "Input/DilutionMock/SangerSequences/")

# ------------------------------------------------------------------------------
# Process Sanger data
# ------------------------------------------------------------------------------

# List all files in the input directory
list.files(path_sanger)

# Store filenames for all 15 species 
Sanger_names <- list.files(
  path_sanger,  pattern = "\\w\\d_S\\d*_", # exclude primer file "27F"
  full.names = TRUE)
Spec_ID <- sub("_.*", "", sub(".*/.._", "", Sanger_names))

# Read files
Sanger_list <- list()
for (i in 1:length(Sanger_names)) {
  Sanger_list[[i]] <- readsangerseq(Sanger_names[i])
}
# Name each Sanger results with the species ID
names(Sanger_list) <- Spec_ID
rm(i, Sanger_names, path_sanger)

# Make base calls (ratio = 0.1 for < 10 CN, 0.05 for 10-20 CN)
for (i in 1:length(Sanger_list)) {
  ratio <- ifelse(names(Sanger_list)[i] %in% c("S02", "S03", "S04"), 0.05, 0.1)
  Sanger_list[[i]] <- makeBaseCalls(Sanger_list[[i]], ratio = ratio)
}
rm(i)

# Order the list
Sanger_list <- Sanger_list[order(names(Sanger_list))]
# Order the ID
Spec_ID <- Spec_ID[order(Spec_ID)]

# Save chromatograms
for (i in 1:length(Sanger_list)) {
  chromatogram(Sanger_list[[i]], width = 100, height = 2, trim5 = 0, trim3 = 0, 
               showcalls = "both", filename = paste0(
                 file_directory, "Output/Dmock_chromatograms/Chroma_", 
                 Spec_ID[i], ".pdf"))
}
rm(i)

# Replace differing basecalls between primary and secondary sequence by .
diff_replace_fun <- function(Sanger_obj) {
  x <- as.character(as.vector(Sanger_obj@primarySeq))
  y <- as.character(as.vector(Sanger_obj@secondarySeq))
  x[which(x != y)] <- "."
  x <- paste0(x, collapse = "")
  return(x)
}
# Apply the function to get the sequences
Sanger_seqs <- unlist(lapply(Sanger_list, diff_replace_fun))

# Count the number of ambiguous basecalls
str_count(Sanger_seqs, "\\.")
str_count(substr(Sanger_seqs[1:15], 68, 430), "\\.")

# Extract the names of expected ASVs, checking in mid part of sequence according to low quality in chromatograms
Expected_seqs <- sapply(substr(Sanger_seqs[1:15], 68, 430), function(x) # 68, 430
  Dmock_ASVtable[grepl(x, Dmock_ASVtable$Sequence), "OTU_ID"])

# Check missing genus level annotation
Expected_seqs
# Missing for ASV_011, ASV_021, ASV_030, which all match "S06" = Escherichia coli

# Extract ASV IDs
Expected_seqs <- unlist(Expected_seqs)
rm(Sanger_list, ratio, Sanger_seqs, Spec_ID, diff_replace_fun)
