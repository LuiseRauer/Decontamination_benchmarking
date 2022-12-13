################################################################################
#
# Decontamination benchmarking on a staggered mock community (Dmock) 
#
################################################################################

# ------------------------------------------------------------------------------
# Set up the workspace
# ------------------------------------------------------------------------------

# Load required packages
library(tidyverse)
# Define directory
file_directory <- "C:/Users/rauerlui/PhD/Projects/02-Decontamination-benchmarking_2020-03/2022.07 Decontamination benchmarking/Decontamination_benchmarking/"

# Define a plot theme
plot_theme <- theme(
  strip.background = element_rect(fill = "#f5f5f5", colour = "grey50"),
  panel.background = element_rect(fill = NA, colour = "grey50"),
  panel.border = element_rect(fill = NA, colour = "grey50"),
  legend.key = element_blank()) 

# Define colours
theme_colours <- c("#4669A5", "#e65d1e", "#5dc9cf", "#ffda0a", "#ba002b", 
                   "#6ed957", "#ffbdc7", "#c97b00", "#fff79c", "#277527", 
                   "#ff7891", "#c6f279", "#3f1163", "#f0b400", "#729ce8", 
                   "#8F1870", "#e3cb6f", "#248c8c", "#e61e1e", "#98d8fe",
                   "#a84c1b", "#D68CCC", "#065e5e", "#ff9c75", "#9DBD32", 
                   "#B957A3", sample(rainbow(250)))

# Define colours for decontamination benchmarking plot
benchm_colours <- 
  c("#EFE11A", "#CDAD0F", "#987C0B", # yellows
    "#FAC771", "#EC921D", "#C45C11", "#6C2C0F", # orange
    "#ECB8FF", "#C583E8", "#7C1CBB", "#2A0140", # purple 
    "#8a8a8a", # grey
    "#A8CE6B", "#50A33B", "#006618", # green
    "#94D8FF", "#37A1DE", "#0061B5", "#00346B") # blue

# ------------------------------------------------------------------------------
# Load and prepare the data
# ------------------------------------------------------------------------------

# Store otu table with count data (samples in rows, ASVs in columns)
otus <- Dmock_ASVtable %>% column_to_rownames("OTU_ID") %>%
  select(-("143-Wasser1":Sequence)) %>%
  t() %>% as.data.frame()
sort(rowSums(otus))
# Remove one sample with only 11 reads
otus <- otus[rowSums(otus) > 1000, ]

# Create OTU table with relative abundances
otus_rel <- sweep(otus, 1, rowSums(otus), "/")
all(round(rowSums(otus_rel), 10) == 1) # Check

# Create metadata
metadata <- read.delim(paste0(
  file_directory, "Input/DilutionMock/Dmock_Metadata.txt"), check.names = FALSE)
# Check that metafile and otu file have the same samples
all(metadata$Sample_ID %in% row.names(otus))
all(row.names(otus) %in% metadata$Sample_ID)

metadata <- data.frame(
  Sample_ID = metadata$Sample_ID,
  Sample_type = ifelse(metadata$Sample_or_Control != "True sample", "NEG2", "SAMPLE"),
  DNA_conc = metadata$quant_reading,
  Input_cells = ifelse(is.na(metadata$Subject), 0, 10^metadata$Subject),
  Dilution = ifelse(metadata$Subject %in% 2:9, paste0(
    "D", chartr("98765432", "01234567", metadata$Subject)), "NEG"))

# Check otu table
sum(otus) # 361651 reads
median(rowSums(otus[row.names(otus) %in% 
                      metadata$Sample_ID[metadata$Sample_type == "SAMPLE"], ]))
median(rowSums(otus[row.names(otus) %in% 
                      metadata$Sample_ID[metadata$Sample_type == "NEG2"], ]))
dim(otus)
length(Expected_seqs)

# Source/sink information for SourceTracker
metadata["SourceSink"] <- ifelse(metadata$Sample_type == "NEG2", "Source", "Sink")

# Load gold standard (classification of ASVs into mock/contaminant by sequence)
Tax_class <- data.frame(
  Contaminant = ifelse(colnames(otus) %in% Expected_seqs, "Mock", "Contaminant"),
  row.names = colnames(otus))

otus_taxa <- Dmock_ASVtable %>% select(Kingdom:Genus)
row.names(otus_taxa) <- Dmock_ASVtable$OTU_ID
# Replace missing or Shigella annotation in sequences that match Escherichia
otus_taxa[row.names(otus_taxa) %in% 
                 c("ASV_011", "ASV_013", "ASV_021", "ASV_030"),
               "Genus"] <- "Escherichia"

# Remove objects from the workspace that are not needed anymore
rm(list = ls()[!ls() %in% c("file_directory", "otus", "otus_rel", "otus_taxa", 
                            "metadata", "Tax_class", "plot_theme", 
                            "theme_colours", "benchm_colours")])

# ------------------------------------------------------------------------------
# Plot sample composition per dilution
# ------------------------------------------------------------------------------

svg(paste0(file_directory, "Output/Plots/Taxonomy_staggered.svg"), 
    width = 5.6, height = 3.8)
otus_rel %>% 
  t() %>% as.data.frame() %>% 
  # Merge with taxonomic information
  merge(., otus_taxa, by = 0, all.x = TRUE) %>% 
  # Separate taxonomy by Gold standard information into mock and contaminant
  mutate(Category = case_when(Row.names %in% row.names(
    Tax_class[Tax_class$Contaminant == "Mock", , drop = FALSE]) ~ Genus,
    TRUE ~ "Contaminants")) %>% 
  # Summarise relative abundances by taxonomy (and sum up all contaminants)
  group_by(Category) %>%
  summarise_at(.vars = metadata$Sample_ID, .funs = sum) %>% 
  # Transform data to long format
  pivot_longer(all_of(metadata$Sample_ID)) %>%
  # Merge with metadata 
  merge(., metadata, by.x = "name", by.y = "Sample_ID", all = TRUE) %>%  # add
  group_by(Category, Dilution) %>% # add
  summarise(value = mean(value)) %>% # add
  ungroup() %>%  # add
  # Add fake row to separate samples and the control by an empty bar
  add_row(Category = "Contaminants", Dilution = "", value = NA) %>%
  mutate(Category = factor(Category, levels = c(unique(sort(.$Category))[c(1:3, 5:12)], "Contaminants")),
         Dilution = factor(Dilution, levels = c(paste0("D", 0:7), "", "NEG"))) %>%
  ggplot(aes(x = Dilution, y = value, fill = Category)) +
  geom_bar(stat = "identity") +
  xlab(bquote('1:10 dilution series ('*10^9*' to'~10^2*' cells)')) + 
  ylab("Relative abundance") +
  plot_theme +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(face = c("plain", rep("bold", 4), rep("plain", 5)),
                                   colour = c("grey30", rep("black", 4), rep("grey30", 5)))) +
  scale_fill_manual("Genus", values = c(
    theme_colours[c(1, 20, 11, 6, 21, 9, 18, 14, 5, 23, 3)], "grey80"))
dev.off()

# ------------------------------------------------------------------------------
# Run SourceTracker
# ------------------------------------------------------------------------------

source(paste0(file_directory, "BM_SourceTracker.R"))
res_sourcetracker_dmock <- BM_sourcetracker(
  file_directory = file_directory, otus = otus, metadata = metadata, 
  Mock_info = Tax_class, raref_depth = 1000)
rm(std.env.colors, eval.fit, jsd, jsdmatrix, kld, plot.eval, 
   plot.sourcetracker.bar, plot.sourcetracker.dist, plot.sourcetracker.fit,
   plot.sourcetracker.pie, predict.sourcetracker, rarefy, run.gibbs, 
   save.mapping.file, sortmatrix, sourcetracker, sourcetracker.error.bars,
   tune.st)

# ------------------------------------------------------------------------------
# Run frequency filter
# ------------------------------------------------------------------------------

source(paste0(file_directory, "BM_Frequency.R"))
res_frequency_dmock <- BM_frequency(
  otus_rel = otus_rel, metadata = metadata, Mock_info = Tax_class, dmax = 7)
rm(BM_frequency)

# ------------------------------------------------------------------------------
# Run NEG2 presence filter
# ------------------------------------------------------------------------------

source(paste0(file_directory, "BM_NegPresence.R"))
res_negpresence_dmock <- BM_negpresence(
  otus_rel = otus_rel, metadata = metadata, Mock_info = Tax_class, dmax = 7)
rm(BM_negpresence)

# ------------------------------------------------------------------------------
# Run MicrobIEM ratio filter
# ------------------------------------------------------------------------------

source(paste0(file_directory, "BM_MicrobIEM.R"))
res_microbiem_dmock <- BM_microbiem(
  otus_rel = otus_rel, metadata = metadata, Mock_info = Tax_class, dmax = 7)
rm(MicrobIEM_decontamination, BM_microbiem)

# ------------------------------------------------------------------------------
# Run Decontam prevalence filter
# ------------------------------------------------------------------------------

source(paste0(file_directory, "BM_DecontamPrev.R"))
res_decontprev_dmock <- BM_decontam_prev(
  otus_rel = otus_rel, metadata = metadata, Mock_info = Tax_class, 
  per_dilution = TRUE, dmax = 7)
rm(BM_decontam_prev)

# ------------------------------------------------------------------------------
# Run Decontam frequency filter
# ------------------------------------------------------------------------------

source(paste0(file_directory, "BM_DecontamFreq.R"))
res_decontfreq_dmock <- BM_decontam_freq(
  otus_rel = otus_rel, metadata = metadata, Mock_info = Tax_class, 
  per_dilution = TRUE, dmax = 7)
rm(BM_decontam_freq)

# ------------------------------------------------------------------------------
# Save the results
# ------------------------------------------------------------------------------

###save(list = ls()[grepl("res_", ls())], file = paste0(file_directory, "Output/R_objects/Dmock_benchm_res.RData"))
###load(paste0(file_directory, "Output/R_objects/Dmock_benchm_res.RData"), verbose = TRUE)

################################################################################
#
# Combination of results 
#
################################################################################

# Combine all benchmarking results
combined_results <- Reduce(
  rbind.data.frame, list(
    data.frame(res_sourcetracker_dmock, Method = "SourceTracker"),
    data.frame(res_frequency_dmock, Method = "Frequency"),
    data.frame(res_negpresence_dmock, Method = "PresenceNEG2"),
    data.frame(res_microbiem_dmock, Method = "MicrobIEM"),
    data.frame(res_decontfreq_dmock, Method = "DecontamFreq"), 
    data.frame(res_decontprev_dmock, Method = "DecontamPrev")))

# Calculate Accuracy, Sensitivity, Specificity, Youden's index, Matthews index
combined_results <- combined_results %>%
  mutate(Accuracy = (TP + TN) / (TP + TN + FP + FN),
         Sensitivity = TP / (TP+FN),
         Specificity = TN / (TN+FP),
         Youden = Sensitivity + Specificity - 1,
         Matthews = (TP*TN - FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
# Correlation between evaluation measures
pairs(combined_results[, c("Accuracy", "Youden", "Matthews")])

# Transform data to long format
combined_results <- combined_results %>% 
  pivot_longer(., cols = c(Accuracy:Matthews), 
               names_to = "Measure", values_to = "Value")

# Add dilution information
combined_results <- merge(combined_results,
                          metadata[, c("Sample_ID", "Dilution")],
                          by.x = "Sample", by.y = "Sample_ID", all.x = TRUE)

# Order methods, filters, and measures for plots
combined_results$Method <- 
  factor(combined_results$Method, levels = c(
    "Frequency", "DecontamFreq", "DecontamPrev", "PresenceNEG2", 
    "SourceTracker", "MicrobIEM"))
combined_results$Measure <- 
  factor(combined_results$Measure, levels = c(
    "Accuracy", "Sensitivity", "Specificity", "Youden", "Matthews"))
combined_results$Filter <- factor(combined_results$Filter, levels = c(
  unique(res_frequency_dmock$Filter),
  unique(res_decontfreq_dmock$Filter),
  unique(res_decontprev_dmock$Filter),
  unique(res_negpresence_dmock$Filter),
  unique(res_sourcetracker_dmock$Filter),
  unique(res_microbiem_dmock$Filter)))

# ------------------------------------------------------------------------------
# Add baseline accuracy (case when not identifying any ASV as contaminant)
# ------------------------------------------------------------------------------

# Define sample IDs
Sample_IDs <- metadata[metadata$Sample_type == "SAMPLE", "Sample_ID"]
# Calculate contaminant prevalence
Cont_prev <- merge(t(otus_rel), Tax_class, by = 0) %>% 
  group_by(Contaminant) %>%
  summarise_at(.vars = Sample_IDs, .funs = sum) %>% 
  filter(Contaminant == "Contaminant") %>% 
  select(all_of(Sample_IDs)) %>% 
  t() %>% as.data.frame() %>%
  mutate(Dilution = paste0(
    "D", chartr("23456789", "76543210", str_sub(row.names(.), -1, -1)))) %>%
  group_by(Dilution) %>%
  summarise(V1 = mean(V1)) 
names(Cont_prev)[2] <- "Cont_prev"
# Add contaminant prevalence to combined results
combined_results <- merge(combined_results, Cont_prev, by = "Dilution", 
                          all.x = TRUE)
# Define 1 as slope for baseline accuracy (when x axis shows contaminant prev.)
combined_results["Cont_prev_line"] <- ifelse(
  combined_results$Measure == "Accuracy", 1, NA)
# Define 0 as slope for "baseline" Youden index
combined_results["Youden_line"] <- ifelse(
  combined_results$Measure == "Youden", 0, NA)
# Define a variable that allows ordering of results for main benchmarking figure
combined_results["X_axis_prox"] <- 
  case_when(grepl("^Frequency", combined_results$Method) ~ 1,
            grepl("^DecontamFreq", combined_results$Method) ~ 2,
            grepl("^DecontamPrev", combined_results$Method) ~ 3,
            grepl("^Prese", combined_results$Method) ~ 4,
            grepl("^Source", combined_results$Method) ~ 5,
            grepl("^MicrobIEM", combined_results$Method) ~ 6)

# ------------------------------------------------------------------------------
# Supp. figure benchmarking
# ------------------------------------------------------------------------------

svg(paste0(file_directory, "Output/Plots/Benchm_staggered_supplement.svg"), 
    width = 9.1, height = 6.5) # width = 12, height = 6
p <- 
  combined_results %>%
  # Calculate mean per dilution (not needed, only for compatibility with Dmock)
  group_by(Filter, Dilution, Measure) %>%
  mutate(Mean = mean(Value),
         Method = factor(Method, levels = c(
           "Frequency", "DecontamFreq", "DecontamPrev", "PresenceNEG2", 
           "SourceTracker", "MicrobIEM"))) %>%
  # Reduce shown filters & select samples
  filter(Dilution %in% c("D1", "D2", "D3", "D4"),
         !Filter %in% c("MicrobIEM, ratio = 1.5", "SourceTracker, beta = 0.001",
                        "Decontam (frequency), thr = 1", 
                        "Decontam (prevalence), thr = 1")) %>% 
  ggplot(aes(x = Cont_prev, y = Mean, colour = Filter, group = Filter)) +
  plot_theme +
  theme(panel.grid.major = element_line(colour = "#e0e0e0"),
        #legend.spacing.y = unit(-0.1, 'cm'),
        legend.position = "none") +
  geom_abline(aes(slope = -Cont_prev_line, intercept = 1), colour = "black") +
  geom_abline(aes(slope = Youden_line, intercept = 0), colour = "black") +
  xlab("Contaminant prevalence") +
  ylab("Mean value per dilution") +
  geom_line(size = 0.7) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(0, 0.8, 0.2),
                     labels = c("0", "0.2", "0.4", "0.6", "0.8")) +
  scale_y_continuous(breaks = seq(0, 1, 0.2),
                     labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0")) +
  scale_colour_manual(values = benchm_colours) +
  facet_grid(Measure ~ Method,
             labeller = labeller(Measure = c("Accuracy" = "Accuracy",
                                             "Sensitivity" = "Sensitivity",
                                             "Specificity" = "Specificity",
                                             "Youden" = "Youden's index",
                                             "Matthews" = "Matthews (MCC)"),
                                 Method = c("Frequency" = "Frequency filter",
                                            "DecontamFreq" = "Decontam (frequency)",
                                            "DecontamPrev" = "Decontam (prevalence)",
                                            "PresenceNEG2" = "Presence filter",
                                            "MicrobIEM" = "MicrobIEM (ratio)",
                                            "SourceTracker" = "SourceTracker"))) +
  guides(colour = guide_legend(ncol = 1, byrow = TRUE)) +
  coord_cartesian(xlim = c(0, 0.8), ylim = c(-0.15, 1.05)) 
p
dev.off()

# ------------------------------------------------------------------------------
# Main figure benchmarking
# ------------------------------------------------------------------------------

svg(paste0(file_directory, "Output/Plots/Benchm_staggered_main.svg"), 
    #width = 12, height = 4.4
    width = 8, height = 3.5)
p <- 
  combined_results %>%
  group_by(Method, Filter, Dilution, Measure, X_axis_prox) %>% 
  summarise(Mean = mean(Value)) %>% 
  filter(Measure %in% c("Youden"),
         !Filter %in% c("MicrobIEM, ratio = 1.5", "SourceTracker, beta = 0.001",
                        "Decontam (frequency), thr = 1", 
                        "Decontam (prevalence), thr = 1"),
         Dilution %in% c("D1", "D2", "D3", "D4")) %>% as.data.frame() %>%
  mutate(Dilution = factor(Dilution, labels = c("D1~(1~x~10^8)",
                                                "D2~(1~x~10^7)",
                                                "D3~(1~x~10^6)",
                                                "D4~(1~x~10^5)"))) %>% 
  ggplot(aes(x = X_axis_prox, y = Mean, colour = Filter, fill = Filter)) +
  plot_theme +
  theme(panel.grid.major.y = element_line(colour = "#e0e0e0"),
        legend.spacing.y = unit(-0.1, 'cm'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(margin = margin(0, 8, 0, 0)),
        legend.position = "none") +
  ylab("Youden's index") +
  geom_vline(xintercept = seq(0.5, 6, by = 1), color = "#e0e0e0", 
             size = 0.5) +
  scale_y_continuous(breaks = seq(0, 1, 0.2),
                     labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0")) +
  #scale_x_discrete(labels = str2expression(paste0("10^{", 9:4, "}"))) +
  geom_col(position = position_dodge2(preserve = "single",
                                      padding = 0.18
                                      )#, width = 0.8
           )+
  scale_colour_manual(name = "Filter \n ", values = benchm_colours) +
  scale_fill_manual(name = "Filter \n ", values = benchm_colours) +
  scale_x_continuous(breaks = seq(1, 6, 1),
                     labels = c("Frequency filter",
                              "Decontam (frequency)",
                              "Decontam (prevalence)",
                              "Presence filter",
                              "SourceTracker",
                              "MicrobIEM (ratio)")) +
  facet_grid(. ~ Dilution, labeller = "label_parsed") +
  coord_cartesian(xlim = c(0.75, 6.25), ylim = c(-0.1, 1)) +
  guides(colour = guide_legend(ncol = 1, byrow = TRUE))
p
dev.off()

# ------------------------------------------------------------------------------
# Text data
# ------------------------------------------------------------------------------

# Contaminant prevalence per sample and NEG
merge(t(otus_rel), Tax_class, by = 0) %>% 
  group_by(Contaminant) %>%
  summarise_at(.vars = metadata$Sample_ID, .funs = sum) %>% 
  filter(Contaminant == "Contaminant") %>% 
  select(all_of(metadata$Sample_ID)) %>% 
  t() %>% as.data.frame()
