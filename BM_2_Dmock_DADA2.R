################################################################################
#
# Microbiome FASTQ denoising with DADA2 (v1.16) for Dmock (staggered mock)
#   based on: https://benjjneb.github.io/dada2/tutorial.html (2021-07-14)
#
################################################################################

# ------------------------------------------------------------------------------
# Set DADA2 parameters
# ------------------------------------------------------------------------------

# Define directory
file_directory <- "C:/Users/rauerlui/PhD/Sonstiges/2022.07 Decontamination benchmarking/Decontamination_benchmarking/"
# Path of input FASTQ files
path_input <- paste0(file_directory, "Input/DilutionMock/FASTQ")
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
# DADA2 installation & packages
# ------------------------------------------------------------------------------

# Install required packages
if (!"dada2" %in% row.names(installed.packages())) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("dada2", version = "3.11")
}

# Install remaining packages
required_packages <- c("tidyverse", "reshape2")
install.packages(setdiff(required_packages, rownames(installed.packages())))
rm(required_packages)

# Load required packages 
library(dada2)
library(tidyverse)
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
filtFs <- file.path(path_output, "Dmock_filtered", 
                    paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path_output, "Dmock_filtered", 
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

# Sort ASVs by decreasing abundance
seqtab.final <- seqtab.final[order(
  rowSums(seqtab.final[, sapply(seqtab.final, is.numeric)]), decreasing = T), ]
rownames(seqtab.final) <- NULL

# Export to Excel/Notepad++/...
write.table(seqtab.final, file = paste0(file_directory, "Input/DilutionMock/Dmock_ASVtable.csv"), 
            quote = FALSE, sep = ";", row.names = FALSE)

# ------------------------------------------------------------------------------
# File export in standard MicrobIEM format
# ------------------------------------------------------------------------------

# Please first run the commands of the previous steps under "Useful format"
# Collapse taxonomy columns
###seqtab.microbIEM <- seqtab.final
###seqtab.microbIEM["Taxonomy"] <- apply(
###  seqtab.final[c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", 
###                 "Species")], 1, function(x) paste(x, collapse = ";"));
###seqtab.microbIEM[c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", 
###                   "Species", "Sequence")] <- NULL
###write.table(seqtab.microbIEM, file = "Proj_ASVtable_MicrobIEM.txt", 
###            quote = FALSE, sep = "\t", row.names = FALSE)
###rm(seqtab.microbIEM)

# ------------------------------------------------------------------------------
# Prepare the data for further analysis with standard scripts
#   (parameter evaluation, merging, metacal)
# ------------------------------------------------------------------------------

# Save point: save all R objects
###save(list = ls(), file = paste0(path_output, "/Dmock_RObjects_DADA2"))

# Tidy up the workspace
rm(list = ls())
