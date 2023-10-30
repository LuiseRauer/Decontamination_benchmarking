################################################################################
#
# Microbiome FASTQ denoising with DADA2 (v1.16) for Staggered mock A &
# ASV classification with Sanger sequencing data
#
################################################################################

# ------------------------------------------------------------------------------
# Set DADA2 parameters
# ------------------------------------------------------------------------------

# Define directory
file_directory <- "C:/Users/rauerlui/PhD/Projects/2020-03_MicrobIEM/2022.07 Decontamination benchmarking/Decontamination_benchmarking/"
# Path of input FASTQ files
path_input <- paste0(file_directory, "Input/Mock_staggered-A_Huelpuesch-Rauer-et-al/FASTQ")
# Path for output
path_output <- paste0(file_directory, "Output")
# Define pattern for distinguishing forward and reverse FASTQ file names
forward_pattern <- "R1_001.fastq"
reverse_pattern <- "R2_001.fastq"

# Set parameters for sequence cutting (V1V3 16S region) and quality filtering
truncLen_par <- c(299, 280)
trimLeft_par <- c(20, 17)
maxEE_par <- c(2, 6) 

# ------------------------------------------------------------------------------
# Set up the workspace
# ------------------------------------------------------------------------------

# Install DADA package
if (!"dada2" %in% row.names(installed.packages())) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("dada2", version = "3.11")
}
# Install sangerSeqR package
if (!"sangerseqR" %in% row.names(installed.packages())) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("sangerseqR")
}

# Install remaining packages
required_packages <- c("tidyverse", "reshape2")
install.packages(setdiff(required_packages, rownames(installed.packages())))
rm(required_packages)

# Load required packages 
library(dada2)
library(tidyverse)
library(sangerseqR)
library(stringr)
library(stringdist)
packageVersion("dada2") # 1.16.0

# Set working directory
setwd(path_output)

# ------------------------------------------------------------------------------
# Load the FASTQ files
# ------------------------------------------------------------------------------

# List all files in the input directory
list.files(path_input)

# Forward and reverse fastq filenames 
fnFs <- sort(list.files(path_input, pattern = forward_pattern, 
                        full.names = TRUE))
fnRs <- sort(list.files(path_input, pattern = reverse_pattern, 
                        full.names = TRUE))
# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# ------------------------------------------------------------------------------
# Inspect read quality profiles
# ------------------------------------------------------------------------------

# Quality profiles for forward reads
plotQualityProfile(fnFs[1:10]) +
  geom_vline(aes(xintercept = trimLeft_par[1]), colour = "blue", alpha = 0.5) +
  geom_vline(aes(xintercept = truncLen_par[1]), colour = "blue", alpha = 0.5)

# Quality profiles for reverse reads
plotQualityProfile(fnRs[1:10]) + 
  geom_vline(aes(xintercept = trimLeft_par[2]), colour = "blue", alpha = 0.5) +
  geom_vline(aes(xintercept = truncLen_par[2]), colour = "blue", alpha = 0.5)

# ------------------------------------------------------------------------------
# Quality filtering and sequence trimming/truncating
# ------------------------------------------------------------------------------

# Create names for filtered files and for a subdirectory of your output folder
filtFs <- file.path(path_output, "Mock_staggered-A_filtered", 
                    paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path_output, "Mock_staggered-A_filtered", 
                    paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Quality filtering
out <- filterAndTrim(
  # R object names:
  fwd = fnFs, filt = filtFs, rev = fnRs, filt.rev = filtRs, 
  # Filtering parameters: 
  truncLen = truncLen_par, trimLeft = trimLeft_par, maxEE = maxEE_par,
  multithread = FALSE, verbose = TRUE)

# ------------------------------------------------------------------------------
# Actual denoising
# ------------------------------------------------------------------------------

# Learn the error rates
errF <- learnErrors(filtFs[file.exists(filtFs)], multithread = FALSE, 
                    verbose = 1)
errR <- learnErrors(filtRs[file.exists(filtRs)], multithread = FALSE, 
                    verbose = 1)

# Visualize the estimated error rates
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

# Sample inference
dadaFs <- dada(filtFs[file.exists(filtFs)], err = errF, multithread = FALSE)
dadaRs <- dada(filtRs[file.exists(filtRs)], err = errR, multithread = FALSE)

# ------------------------------------------------------------------------------
# Merging & chimera removal
# ------------------------------------------------------------------------------

# Merge forward and reverse reads
mergers <- mergePairs(dadaFs, filtFs[file.exists(filtFs)], dadaRs, 
                      filtRs[file.exists(filtRs)], verbose = TRUE)

# Construct a sequence table for all merged reads
seqtab <- makeSequenceTable(mergers)

# Chimera removal
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", 
                                    multithread = TRUE, verbose = TRUE)

# ------------------------------------------------------------------------------
# Sanity checks
# ------------------------------------------------------------------------------

# Number of samples & Number of ASVs
dim(seqtab.nochim) 
# Distribution of sequence lengths
hist(nchar(getSequences(seqtab.nochim)))
# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- list(data.frame(out, row.names = sample.names), 
              data.frame(sapply(dadaFs, getN)), 
              data.frame(sapply(dadaRs, getN)), 
              data.frame(sapply(mergers, getN)), 
              data.frame(rowSums(seqtab.nochim)))
track <- lapply(track, function(x) rownames_to_column(x, "sample"))
track <- track %>% reduce(merge, by = "sample", all = TRUE)
track <- track %>% column_to_rownames(var = "sample")
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
# Mean percentage of recovered reads per step
round(colMeans(track/track$input*100, na.rm = TRUE), 1)

################################################################################
# 
# Additional steps
#
################################################################################

# ------------------------------------------------------------------------------
# Additional quality control
# ------------------------------------------------------------------------------

# Create an ASV table for filtering steps
seqtab.filt <- seqtab.nochim

# Filter sequences by length (+/- 20% of expected length, i.e. 526 bp)
seqtab.filt <- seqtab.filt[, nchar(colnames(seqtab.filt)) %in% 421:631]

# Singleton removal not needed
table(apply(seqtab.filt, 2, max) > 1) # 293 TRUE

# Create an ID for each ASV
ASV_ID <- paste0("ASV_", formatC(seq(ncol(seqtab.filt)), 
                                 width = nchar(ncol(seqtab.filt)), flag = "0"))

# ------------------------------------------------------------------------------
# Taxonomic annotation through DADA2 with RDP
# ------------------------------------------------------------------------------

# Download the file RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz from 
#   https://zenodo.org/record/4735821#.YO9BsUHgo2w and save the database in 
#   your output directory
Sys.time() # takes 15min
set.seed(1)
taxa_RDP <- assignTaxonomy(
  seqtab.filt, 
  paste0(file_directory, "Input/RDPClassifier/RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz"),
  multithread = FALSE)
Sys.time()

# Merge ASV table with RDP taxonomy
seqtab.final <- merge(data.frame(OTU_ID = ASV_ID, t(seqtab.filt),
                                 check.names = FALSE), 
                      taxa_RDP, by = 0, all = TRUE)

# ------------------------------------------------------------------------------
# Useful format for export to Excel or to continue work in R
# ------------------------------------------------------------------------------

# Move the sequence column to the end
seqtab.final["Sequence"] <- seqtab.final["Row.names"]
# Remove merging column
seqtab.final["Row.names"] <- NULL

# Sort ASVs by OTU_ID
seqtab.final <- seqtab.final[order(seqtab.final$OTU_ID), ]
rownames(seqtab.final) <- NULL

# Export to Excel
write.table(seqtab.final, file = paste0(
  file_directory, "Input/Mock_staggered-A_Huelpuesch-Rauer-et-al/Mock_staggered-A_ASVs.csv"), 
  quote = FALSE, sep = ";", row.names = FALSE)

# ------------------------------------------------------------------------------
# Save results
# ------------------------------------------------------------------------------

# Save point: save all R objects
###save(list = ls(), file = paste0(path_output, "/R_objects/0_DADA2_staggered-A.RData"))
###load(paste0(path_output, "/R_objects/0_DADA2_staggered-A.RData"), verbose = TRUE)

################################################################################
#
# Classify ASVs with Sanger sequencing data 
#
################################################################################

# ------------------------------------------------------------------------------
# Process Sanger data
# ------------------------------------------------------------------------------

# Path of input FASTQ files (must be demultiplexed, can be fastq.gz)
path_sanger <- paste0(file_directory, "Input/Mock_staggered-A_Huelpuesch-Rauer-et-al/SangerSequences/")

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
                 file_directory, "Output/Mock_staggered-A_chromatograms/Chroma_", 
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
# Count the number of ambiguous basecalls in the high-quality middle region
str_count(substr(Sanger_seqs, 68, 430), "\\.")

# Extract expected ASVs that match the high-quality mid part of Sanger sequences
Expected_seqs_list <- sapply(substr(Sanger_seqs, 68, 430), function(x) # 68, 430
  seqtab.final[grepl(x, seqtab.final$Sequence), "OTU_ID"])

# Check expected ASVs (S15 (S. pneumoniae) not detected in the mock)
Expected_seqs_list

# Unlist ASV IDs
Expected_seqs <- unlist(Expected_seqs_list)

# ------------------------------------------------------------------------------
# Add sequences with small sequence errors to expected sequences
# ------------------------------------------------------------------------------

# Calculate LV distance between true sequences and remaining sequences
Exp_seq_dist <- stringdistmatrix(
  # Remaining ASV sequences
  setNames(seqtab.final[!seqtab.final$OTU_ID %in% Expected_seqs, "Sequence"],
           seqtab.final[!seqtab.final$OTU_ID %in% Expected_seqs, "OTU_ID"]),
  # Sequences of expected ASVs
  setNames(seqtab.final[seqtab.final$OTU_ID %in% Expected_seqs, "Sequence"],
           seqtab.final[seqtab.final$OTU_ID %in% Expected_seqs, "OTU_ID"]),
  useNames = "names")
# Extract the minimum LV distance per row (remaining ASV sequences)
Min_seq_dist <- sort(apply(Exp_seq_dist, 1, min))
# Add ASVs with <= 4 LV distance to other expected ASVs
Expected_seqs <- c(Expected_seqs, names(Min_seq_dist[Min_seq_dist <= 4]))
# Is RDP taxonomic annotation correct for additional sequences? Yes.
View(as.matrix(Exp_seq_dist)[row.names(as.data.frame(Exp_seq_dist)) %in% 
                    names(Min_seq_dist[Min_seq_dist <= 4]), ])
Expected_seqs_list
# ASV_252: S04, ASV_167: S01, ASV_143: S03, ASV_118: S03, ASV_111: S09, 
#   ASV_109: S06, ASV_100: S13
seqtab.final[seqtab.final$OTU_ID %in% Expected_seqs & 
               !seqtab.final$OTU_ID %in% unlist(Expected_seqs_list), 
             c("OTU_ID", "Genus")]

# ------------------------------------------------------------------------------
# Check lost S. pneumoniae
# ------------------------------------------------------------------------------

# Distance between expected sequences of S. pneumoniae and S. mutans
stringdist(substr(Sanger_seqs, 68, 430)[14], substr(Sanger_seqs, 68, 430)[15])

# Calculate LV distance between true S. pneum. sequence and remaining sequences
Exp_seq_dist_Spneu <- stringdistmatrix(
  # Remaining ASV sequences, cut to the the high-quality Sanger sequ. subregion
  setNames(substr(
    seqtab.final[!seqtab.final$OTU_ID %in% Expected_seqs, "Sequence"], 54, 416),
    seqtab.final[!seqtab.final$OTU_ID %in% Expected_seqs, "OTU_ID"]),
  # Sequence of S. pneumoniae
  setNames(substr(Sanger_seqs[15], 68, 430), "S15"),
  useNames = "names")
# Extract the minimum LV distance per row (remaining ASV sequences)
sort(apply(Exp_seq_dist_Spneu, 1, min)) # all >= 50

# ------------------------------------------------------------------------------
# Check missing genus level annotation
# ------------------------------------------------------------------------------

Missing_annot <- seqtab.final[seqtab.final$OTU_ID %in% Expected_seqs & 
                                is.na(seqtab.final$Genus), "OTU_ID"]
lapply(Expected_seqs_list, function(x) Missing_annot %in% x) 
# All match S06 = E. coli

# Is RDP taxonomic annotation correct for expected sequences? Yes.
View(cbind.data.frame(
  seqtab.final[seqtab.final$OTU_ID %in% Expected_seqs, "OTU_ID"],
  seqtab.final[seqtab.final$OTU_ID %in% Expected_seqs, "Genus"]))
Expected_seqs_list

# ------------------------------------------------------------------------------
# Export files needed for benchmarking
# ------------------------------------------------------------------------------

# Create metadata
metadata <- read.delim(paste0(
  file_directory, "Input/Mock_staggered-A_Huelpuesch-Rauer-et-al/Mock_staggered-A_Metadata.txt"), check.names = FALSE)

# Save point: save all R objects
###save(seqtab.final, Expected_seqs, metadata, file = paste0(file_directory, "/Input/Mock_staggered-A_Huelpuesch-Rauer-et-al/Mock_staggered-A.RData"))
###load(paste0(file_directory, "/Input/Mock_staggered-A_Huelpuesch-Rauer-et-al/Mock_staggered-A.RData"), verbose = TRUE)

# Tidy up the workspace
rm(list = ls())
