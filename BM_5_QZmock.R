################################################################################
#
# Decontamination benchmarking of the QZmock
#
################################################################################

# ------------------------------------------------------------------------------
# Set up the workspace
# ------------------------------------------------------------------------------

# Load required packages
library(tidyverse)
library(cowplot) # plot_grid, get_legend
library(ggpubr) # as_ggplot

# Define directory
file_directory <- "C:/Users/rauerlui/PhD/Projects/2020-03_MicrobIEM/2022.07 Decontamination benchmarking/Decontamination_benchmarking/"

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
                   "#B957A3", "#83D9AB", "#D32868", "#FFCC85", "#B9A8F0",
                   sample(rainbow(250)))

# Define colours for decontamination benchmarking plot
benchm_colours <- 
  c("#EFE11A", "#CDAD0F", "#987C0B", # yellow - frequency
    "#FAC771", "#EC921D", "#C45C11", "#6C2C0F", # orange - Decontam freq
    "#ECB8FF", "#C583E8", "#7C1CBB", "#2A0140", # purple - Decontam prev
    "#A8CE6B", "#50A33B", "#006618", # green - SourceTracker
    "#333333", ##8a8a8a", # grey - presence filter
    "#C0C0C0", "#8E8E8E", "#5C5C5C", "#333333", # grey - MicrobIEM span
    "#94D8FF", "#37A1DE", "#0061B5", "#00346B") # blue - MicrobIEM ratio


benchm_comb_colours <- c("#A9C8DA", "#6E9EB9", "#3F6383", "#183553", "black")

# ------------------------------------------------------------------------------
# Load and prepare the data
# ------------------------------------------------------------------------------

# Define otu table with count data (samples in rows, ASVs in columns)
load(paste0(file_directory, "Input/Mock_spike_Rauer-de-Tomassi-et-al/SpikeMock_ASVtable.RData"), 
     verbose = TRUE)
otus_orig <- otus
otus <- otus %>% 
  column_to_rownames("ASV_ID") %>% 
  select(-LV_4) %>%
  t(.) %>% as.data.frame(.)

# Create OTU table with relative abundances
otus_rel <- sweep(otus, 1, rowSums(otus), "/")
all(rowSums(otus_rel) == 1) # Check

# Create metadata
load(paste0(file_directory, "Input/Mock_spike_Rauer-de-Tomassi-et-al/SpikeMock_Metadata.RData"), 
     verbose = TRUE)
# Check that metafile and otu file have the same samples
all(metadata$Sample_ID %in% row.names(otus))
all(row.names(otus) %in% metadata$Sample_ID)

# Add source/sink information for SourceTracker
metadata["SourceSink"] <- ifelse(metadata$Sample_type == "NEG2", "Source", 
                                 "Sink")

# Create gold standard (classify ASVs into mock/contaminant by sequence)
mock_taxa <- otus_orig[grepl("Truep|Allob|Imtech", otus_orig$LV_4), "ASV_ID"]
# Add cross-contamination from skin:
cross_cont_taxa <- c(otus_orig[grepl("_0003|_0045", otus_orig$ASV_ID), "ASV_ID"])
# Add cross-contamination from other mock samples:
cross_cont_taxa <- c(cross_cont_taxa, otus_orig[!is.na(otus_orig$LV_4), "ASV_ID"])

Tax_class <- data.frame(
  Contaminant = ifelse(colnames(otus) %in% mock_taxa, "Mock", "Contaminant"),
  row.names = colnames(otus))

# Store taxonomic information
otus_taxa <- otus_orig %>% select(LV_4)
row.names(otus_taxa) <- otus_orig$ASV_ID
otus_taxa["Genus"] <- gsub("_.*", "", otus_taxa$LV_4, perl = TRUE)

# ------------------------------------------------------------------------------
# Plot sample composition per dilution
# ------------------------------------------------------------------------------

svg(paste0(file_directory, "Output/Plots/F02_Taxonomy_spike.svg"), 
    width = 3.75, height = 3.8)
otus_rel %>% 
  t() %>% as.data.frame() %>% 
  # Merge with taxonomic information
  merge(., otus_taxa, by = 0, all.x = TRUE) %>% 
  # Separate taxonomy by Gold standard information into mock and contaminant
  mutate(Category = case_when(
    grepl("Allob|Imtech|Truep", LV_4) & !is.na(Genus) ~ Genus,
    Row.names %in% cross_cont_taxa ~ "Cross-contaminants",
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
  mutate(Category = factor(Category, levels = c(unique(sort(.$Category))[c(1, 4, 5)], "Cross-contaminants", "Contaminants")),
         Dilution = factor(Dilution, levels = c(paste0("D", 0:1), "", "NEG"))) %>%
  ggplot(aes(x = Dilution, y = value, fill = Category)) +
  geom_bar(stat = "identity") +
  xlab(bquote('1:20 dilution ('*10^5*' to 6*'*10^3*' cells)')) + 
  ylab("Relative abundance") + ggtitle("Staggered mock community B") +
  plot_theme +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.text.align = 0,
        plot.title = element_text(size = rel(1.1), hjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(face = c(rep("bold", 2), rep("plain", 2)),
                                   colour = c(rep("black", 2), rep("grey30", 2)))) +
  scale_fill_manual("Genus", values = c(theme_colours[27:29], "grey50", "grey80"),
    labels = c("Contaminants" = "Contaminants",
               "Cross-contaminants" = "Cross-contaminants",
               "Allobacillus" = expression(italic("Allobacillus")),
               "Imtechella" = expression(italic("Imtechella")),
               "Truepera" = expression(italic("Truepera"))))
dev.off()

# ------------------------------------------------------------------------------
# Run SourceTracker
# ------------------------------------------------------------------------------

source(paste0(file_directory, "BM_SourceTracker.R"))
res_sourcetracker_spike <- BM_sourcetracker(
  file_directory = file_directory, otus = otus, metadata = metadata, 
  Mock_info = Tax_class, raref_depth = 1000) # change me later to 10000
rm(std.env.colors, eval.fit, jsd, jsdmatrix, kld, plot.eval, 
   plot.sourcetracker.bar, plot.sourcetracker.dist, plot.sourcetracker.fit,
   plot.sourcetracker.pie, predict.sourcetracker, rarefy, run.gibbs, 
   save.mapping.file, sortmatrix, sourcetracker, sourcetracker.error.bars,
   tune.st, BM_sourcetracker)

# ------------------------------------------------------------------------------
# Run frequency filter
# ------------------------------------------------------------------------------

source(paste0(file_directory, "BM_Frequency.R"))
res_frequency_spike <- BM_frequency(
  otus_rel = otus_rel, metadata = metadata, Mock_info = Tax_class, dmax = 1)
rm(BM_frequency)

# ------------------------------------------------------------------------------
# Run NEG2 presence filter
# ------------------------------------------------------------------------------

source(paste0(file_directory, "BM_NegPresence.R"))
res_negpresence_spike <- BM_negpresence(
  otus_rel = otus_rel, metadata = metadata, Mock_info = Tax_class, dmax = 1)
rm(BM_negpresence)

# ------------------------------------------------------------------------------
# Run MicrobIEM ratio filter
# ------------------------------------------------------------------------------

source(paste0(file_directory, "BM_MicrobIEM_ratio.R"))
res_microbiemratio_spike <- BM_microbiem_ratio(
  otus_rel = otus_rel, metadata = metadata, Mock_info = Tax_class, dmax = 1)
rm(MicrobIEM_decontamination, BM_microbiem_ratio)

# ------------------------------------------------------------------------------
# Run MicrobIEM span filter
# ------------------------------------------------------------------------------

source(paste0(file_directory, "BM_MicrobIEM_span.R"))
res_microbiemspan_spike <- BM_microbiem_span(
  otus_rel = otus_rel, metadata = metadata, Mock_info = Tax_class, dmax = 1)
rm(MicrobIEM_decontamination, BM_microbiem_span)

# ------------------------------------------------------------------------------
# Run Decontam prevalence filter
# ------------------------------------------------------------------------------

source(paste0(file_directory, "BM_Decontam_prev.R"))
res_decontprev_spike <- BM_decontam_prev(
  otus_rel = otus_rel, metadata = metadata, Mock_info = Tax_class, 
  per_dilution = TRUE, dmax = 1)
rm(BM_decontam_prev)

# ------------------------------------------------------------------------------
# Run Decontam frequency filter
# ------------------------------------------------------------------------------

source(paste0(file_directory, "BM_Decontam_freq.R"))
res_decontfreq_spike <- BM_decontam_freq(
  otus_rel = otus_rel, metadata = metadata, Mock_info = Tax_class, 
  per_dilution = TRUE, dmax = 1)
rm(BM_decontam_freq)

# ------------------------------------------------------------------------------
# Save the results
# ------------------------------------------------------------------------------

###save(list = ls()[grepl("res_", ls())], file = paste0(file_directory, "Output/R_objects/3_BM_spike_res.RData"))
###load(paste0(file_directory, "Output/R_objects/3_BM_spike_res.RData"), verbose = TRUE)

################################################################################
#
# Combination of results 
#
################################################################################

# Combine all benchmarking results
combined_results <- Reduce(
  rbind.data.frame, list(
    data.frame(res_sourcetracker_spike, Method = "SourceTracker"),
    data.frame(res_frequency_spike, Method = "Frequency"),
    data.frame(res_negpresence_spike, Method = "PresenceNEG2"),
    data.frame(res_microbiemratio_spike, Method = "MicrobIEMRatio"),
    data.frame(res_microbiemspan_spike, Method = "MicrobIEMSpan"),
    data.frame(res_decontfreq_spike, Method = "DecontamFreq"), 
    data.frame(res_decontprev_spike, Method = "DecontamPrev")))

# Calculate Accuracy, Sensitivity, Specificity, Youden's index, Matthews index
combined_results <- combined_results %>%
  mutate(Accuracy = (TP + TN) / (TP + TN + FP + FN),
         Sensitivity = TP / (TP+FN),
         Specificity = TN / (TN+FP),
         Precision = TP / (TP+FP),
         Youden = Sensitivity + Specificity - 1,
         Matthews = (TP*TN - FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
# Replace NaN values in MCC
combined_results <- combined_results %>%
  mutate(Matthews = case_when(
    (is.na(Matthews) & TP == 0 & FP == 0 & FN > 0 & TN > 0) ~ 0,
    TRUE ~ Matthews),
    Precision = case_when(
      (is.na(Precision) & TP == 0 & (FP > 0 | FN > 0)) ~ 0,
      TRUE ~ Precision)
  )

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
    "Frequency", "DecontamFreq", "DecontamPrev", "SourceTracker", 
    "PresenceNEG2", "MicrobIEMSpan", "MicrobIEMRatio"))
combined_results$Measure <- 
  factor(combined_results$Measure, levels = c(
    "Accuracy", "Sensitivity", "Specificity", "Precision", "Youden", "Matthews"))
combined_results$Filter <- factor(combined_results$Filter, levels = c(
  unique(res_frequency_spike$Filter),
  unique(res_decontfreq_spike$Filter),
  unique(res_decontprev_spike$Filter),
  unique(res_sourcetracker_spike$Filter),
  unique(res_negpresence_spike$Filter),
  #unique(res_microbiemspan_spike$Filter),
  rev(c("MicrobIEM, span = 0.25", "MicrobIEM, span = 0.5", "MicrobIEM, span = 0.75",
    "MicrobIEM, span = 1")),
  unique(res_microbiemratio_spike$Filter)))

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
  merge(., metadata[, c("Sample_ID", "Dilution")], by.x = 0, by.y = "Sample_ID", 
        all.x = TRUE) %>% #View
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
# Define 0 as slope for "baseline" MCC
combined_results["Matthews_line"] <- ifelse(
  combined_results$Measure == "Matthews", 0, NA)

# Define a variable that allows ordering of results for main benchmarking figure
combined_results["X_axis_prox"] <- 
  case_when(grepl("^Frequency", combined_results$Method) ~ 1,
            grepl("^DecontamFreq", combined_results$Method) ~ 2,
            grepl("^DecontamPrev", combined_results$Method) ~ 3,
            grepl("^Source", combined_results$Method) ~ 4,
            grepl("^Prese", combined_results$Method) ~ 5,
            grepl("^MicrobIEMSpan", combined_results$Method) ~ 6,
            grepl("^MicrobIEMRatio", combined_results$Method) ~ 7)

# ------------------------------------------------------------------------------
# Supp. figure benchmarking
# ------------------------------------------------------------------------------

svg(paste0(file_directory, "Output/Plots/S04_Benchm_spike_supp.svg"), 
    width = 9, height = 7.9) # width = 12, height = 6
p <- 
  combined_results %>%
  # Calculate mean per dilution (not needed, only for compatibility with Dmock)
  group_by(Filter, Dilution, Measure) %>%
  mutate(Mean = mean(Value)) %>%
  # Reduce shown filters & select samples
  filter(!Filter %in% c("MicrobIEM, ratio = 1.5", "SourceTracker, beta = 0.001",
                        "Decontam (frequency), thr = 1", 
                        "Decontam (prevalence), thr = 1")) %>% 
  ggplot(aes(x = Cont_prev, y = Mean, colour = Filter, group = Filter)) +
  plot_theme +
  theme(panel.grid.major = element_line(colour = "#e0e0e0"),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(0, 0, 0, 0)), 
        plot.title = element_text(size = rel(1.1))) +
  geom_abline(aes(slope = -Cont_prev_line, intercept = 1), colour = "black") +
  geom_abline(aes(slope = Youden_line, intercept = 0), colour = "black") +
  geom_abline(aes(slope = Matthews_line, intercept = 0), colour = "black") +
  xlab("Contaminant prevalence") +
  ylab("Mean value per dilution") + ggtitle("Staggered mock community B") +
  geom_line(size = 0.7) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(0, 0.8, 0.2),
                     labels = c("0", "0.2", "0.4", "0.6", "0.8")) +
  scale_y_continuous(breaks = seq(-0.2, 1, 0.2),
                     labels = c("-0.2", "0", "0.2", "0.4", "0.6", "0.8", "1.0")) +
  scale_colour_manual(values = benchm_colours) +
  facet_grid(Measure ~ Method,
             labeller = labeller(Measure = c("Accuracy" = "Accuracy",
                                             "Sensitivity" = "Sensitivity",
                                             "Specificity" = "Specificity",
                                             "Precision" = "Precision",
                                             "Youden" = "Youden's index",
                                             "Matthews" = "Matthews (MCC)"),
                                 Method = c("Frequency" = "Frequency filter",
                                            "DecontamFreq" = "Decontam (freq.)",
                                            "DecontamPrev" = "Decontam (prev.)",
                                            "PresenceNEG2" = "Presence filter",
                                            "MicrobIEMRatio" = "MicrobIEM (ratio)",
                                            "MicrobIEMSpan" = "MicrobIEM (span)",
                                            "SourceTracker" = "SourceTracker"))) +
  guides(colour = guide_legend(ncol = 1, byrow = TRUE)) +
  coord_cartesian(xlim = c(0, 0.8), ylim = c(-0.3, 1.05)) 
p
dev.off()

# ------------------------------------------------------------------------------
# Main figure benchmarking
# ------------------------------------------------------------------------------

svg(paste0(file_directory, "Output/Plots/F03_Benchm_spike_main.svg"), 
    #width = 12, height = 4.4
    width = 4.29, height = 3.2)
p <- 
  combined_results %>%
  group_by(Method, Filter, Dilution, Measure, X_axis_prox) %>% 
  summarise(Mean = mean(Value)) %>% 
  filter(Measure %in% c("Youden"),
         !Filter %in% c("MicrobIEM, ratio = 1.5", "SourceTracker, beta = 0.001",
                        "Decontam (frequency), thr = 1", 
                        "Decontam (prevalence), thr = 1"),
         #Dilution %in% c("D1", "D2", "D3", "D4")
         ) %>% 
  as.data.frame() %>%
  mutate(Dilution = factor(Dilution, labels = c("D0~(1~x~10^5)",
                                                "D1~(6~x~10^3)"))) %>% 
  ggplot(aes(x = X_axis_prox, y = Mean, colour = Filter, fill = Filter)) +
  plot_theme +
  theme(panel.grid.major.y = element_line(colour = "#e0e0e0"),
        legend.spacing.y = unit(-0.1, 'cm'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = rel(1.1)),
        axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
        legend.position = "none"
        ) +
  ylab("Youden's index") + ggtitle("Staggered mock community B") +
  geom_vline(xintercept = seq(0.5, 7, by = 1), color = "#e0e0e0", 
             size = 0.5) +
  scale_y_continuous(breaks = seq(-0.2, 1, 0.2),
                     labels = c("-0.2", "0", "0.2", "0.4", "0.6", "0.8", "1.0")) +
  #scale_x_discrete(labels = str2expression(paste0("10^{", 9:4, "}"))) +
  geom_col(position = position_dodge2(preserve = "single",
                                      padding = 0.18))+
  scale_colour_manual(name = "Filter \n ", values = benchm_colours) +
  scale_fill_manual(name = "Filter \n ", values = benchm_colours) +
  scale_x_continuous(breaks = seq(1, 7, 1),
                     labels = c("Frequency filter",
                                "Decontam (freq.)",
                                "Decontam (prev.)",
                                "SourceTracker",
                                "Presence filter",
                                "MicrobIEM (span)",
                                "MicrobIEM (ratio)")) +
  facet_grid(. ~ Dilution, labeller = "label_parsed") +
  coord_cartesian(xlim = c(0.75, 7.25), ylim = c(-0.3, 1)) +
  guides(colour = guide_legend(ncol = 1, byrow = TRUE))
p
dev.off()

# ------------------------------------------------------------------------------
# Legend
# ------------------------------------------------------------------------------

# Build legend per filter tool
for (i in c("MicrobIEMRatio", "MicrobIEMSpan", "PresenceNEG2", "SourceTracker", 
            "DecontamPrev", "DecontamFreq", "Frequency"
)) {
  p <- 
    combined_results %>%
    filter(Method == i) %>% 
    filter(!Filter %in% c("MicrobIEM, ratio = 1.5", "SourceTracker, beta = 0.001",
                          "Decontam (frequency), thr = 1", 
                          "Decontam (prevalence), thr = 1"),
           Dilution %in% c("D0", "D1")
    ) %>% as.data.frame() %>%
    # Replace tool name in legend
    mutate(Filter = gsub(".*, ", "", Filter)) %>% 
    # Prevent alphabetic ordering in MicrobIEM ratio and SourceTracker
    {if (i %in% c("MicrobIEMRatio", "SourceTracker")) 
      mutate(., Filter = factor(Filter, levels = gsub(
        ".*, ", "", levels(combined_results$Filter))[c(14:16, 23, 25:27)])) else .} %>%
    ggplot(aes(x = X_axis_prox, y = Value, colour = Filter, fill = Filter)) +
    plot_theme +
    theme(legend.spacing.y = unit(-0.1, 'cm'),
          legend.position = "left",
          legend.title = element_text(margin = margin(0, 0, 6, 0)),
          legend.justification = c(0, 0.8)) +
    geom_col(position = position_dodge2(preserve = "single")) +
    guides(colour = guide_legend(ncol = 1, byrow = TRUE))
  if (i == "MicrobIEMRatio") {
    p1 <- p + scale_colour_manual(name = "MicrobIEM (ratio)", values = benchm_colours[20:23]) +
      scale_fill_manual(name = "MicrobIEM (ratio)", values = benchm_colours[20:23])
  }
  if (i == "MicrobIEMSpan") {
    p2 <- p + scale_colour_manual(name = "MicrobIEM (span)", values = benchm_colours[16:19],
                                  labels = c("span = 4 of all", "span = 3 of all", 
                                             "span = 2 of all", "span = 1 of all")
                                  ) +
      scale_fill_manual(name = "MicrobIEM (span)", values = benchm_colours[16:19],
                        labels = c("span = 4 of all", "span = 3 of all", 
                                   "span = 2 of all", "span = 1 of all")
                        )
  }
  if (i == "PresenceNEG2") {
    p3 <- p + scale_colour_manual(name = "Presence filter", values = benchm_colours[15]) +
      scale_fill_manual(name = "Presence filter", values = benchm_colours[15])
  }
  if (i == "SourceTracker") {
    p4 <- p + scale_colour_manual(name = "SourceTracker", values = benchm_colours[12:14]) +
      scale_fill_manual(name = "SourceTracker", values = benchm_colours[12:14])
  }
  if (i == "DecontamPrev") {
    p5 <- p + scale_colour_manual(name = "Decontam (prev.)", values = benchm_colours[8:11]) +
      scale_fill_manual(name = "Decontam (prev.)", values = benchm_colours[8:11])
  }
  if (i == "DecontamFreq") {
    p6 <- p + scale_colour_manual(name = "Decontam (freq.)", values = benchm_colours[4:7]) +
      scale_fill_manual(name = "Decontam (freq.)", values = benchm_colours[4:7])
  }
  if (i == "Frequency") {
    p7 <- p + scale_colour_manual(name = "Frequency filter", values = benchm_colours) +
      scale_fill_manual(name = "Frequency filter", values = benchm_colours)
  }
}

svg(paste0(file_directory, "Output/Plots/Benchm_legend.svg"),
    width = 1.8, height = 7)
plot_grid(as_ggplot(get_legend(p7)), as_ggplot(get_legend(p6)), 
          as_ggplot(get_legend(p5)), as_ggplot(get_legend(p4)), 
          as_ggplot(get_legend(p3)), as_ggplot(get_legend(p2)), 
          as_ggplot(get_legend(p1)), ncol = 1, 
          align = "hv", rel_heights = c(0.91, 1.1, 1.1, 0.91, 0.55, 1.1, 1.1))
dev.off()

# ------------------------------------------------------------------------------
# Benchmarking of combination of MicrobIEM ratio and span filter
# ------------------------------------------------------------------------------

# Check the amount of cross-contamination into negative controls
merge(t(otus_rel), Tax_class, by = 0, all = TRUE) %>% 
  filter(Contaminant == "Mock") %>%
  View 
# 2/4 negative controls contains 1 and 2 mock ASVs

# Run combination of MicrobIEM ratio and span filter
source(paste0(file_directory, "BM_MicrobIEM_comb.R"))
res_microbiem_comb_spike <- BM_microbiem_comb(
  otus_rel = otus_rel, metadata = metadata, Mock_info = Tax_class, dmax = 1)

# Add dilution information
res_microbiem_comb <- merge(
  res_microbiem_comb_spike,
  metadata[, c("Sample_ID", "Dilution")],
  by.x = "Sample", by.y = "Sample_ID", all.x = TRUE)

svg(paste0(file_directory, "Output/Plots/S03_Benchm_spike_MicrobIEM_supp.svg"), 
    #width = 12, height = 4.4
    width = 6.49, height = 4.05)
res_microbiem_comb %>% 
  #filter(Sample != "550-BiomPsy-979") %>%
  separate(Filter, c("Ratio", "Span"), "; ") %>% 
  filter(!Ratio %in% c("MicrobIEM, ratio = 1.5")) %>%
  mutate(Accuracy = (TP + TN) / (TP + TN + FP + FN),
         Sensitivity = TP / (TP+FN),
         Specificity = TN / (TN+FP),
         Youden = Sensitivity + Specificity - 1,
         Matthews = (TP*TN - FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)),
         Span2 = recode(Span, "span = NA" = "span = 0 of all (no span filter)", 
                        "span = 0.25" = "span = 1 of all",
                        "span = 0.5" = "span = 2 of all",
                        "span = 0.75" = "span = 3 of all",
                        "span = 1" = "span = 4 of all"),
         Span2 = factor(Span2, levels = rev(c("span = 0 of all (no span filter)", 
                                              "span = 1 of all", "span = 2 of all", 
                                              "span = 3 of all", "span = 4 of all"))),
         Ratio = gsub("MicrobIEM, ", "", Ratio),
         Ratio = recode(Ratio, "ratio = NA" = "ratio = NA (no ratio filter)",),
         Ratio = factor(Ratio, levels = c("ratio = NA (no ratio filter)",
                                          "ratio = 0.1", "ratio = 0.5",
                                          "ratio = 1", "ratio = 2"))) %>% 
  mutate(Dilution = factor(Dilution, labels = c("D0~(1~x~10^5)",
                                                "D1~(6~x~10^3)"))) %>% 
  group_by(Dilution, Span2, Ratio) %>%
  summarise(Youden = mean(Youden)) %>%
  #View
  ggplot(., aes(x = Ratio, y = Youden, fill = Span2, colour = Span2)) +
  geom_col(position = position_dodge2(preserve = "single",
                                      padding = 0.18))+
  facet_grid(. ~ Dilution, labeller = "label_parsed") +
  xlab("MicrobIEM (ratio)") + ylab("Youden's index") + 
  ggtitle("Staggered mock community B") +
  scale_y_continuous(breaks = seq(0, 1, 0.2),
                     labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0")) +
  coord_cartesian(ylim = c(-0.1, 1)) +
  scale_fill_manual("MicrobIEM (span)", values = benchm_comb_colours) +
  scale_colour_manual("MicrobIEM (span)", values = benchm_comb_colours) +
  plot_theme +
  theme(panel.grid.major.y = element_line(colour = "#e0e0e0"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = rel(1.1)),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(margin = margin(0, 8, 0, 0))) +
  geom_vline(xintercept = seq(0.5, 6, by = 1), color = "#e0e0e0", 
             size = 0.5)
dev.off()

####
library(tidyverse)
# Contaminant prevalence per sample and NEG
merge(t(otus_rel), Tax_class, by = 0) %>% 
  group_by(Contaminant) %>% #View
  summarise_at(.vars = c(metadata$Sample_ID), .funs = sum) %>% 
  filter(Contaminant == "Contaminant") %>% 
  select(all_of(metadata$Sample_ID)) %>% 
  t() %>% as.data.frame() %>% arrange(V1)
# 1.8-55.5 % contamination

median(c(0.22969236, 0.25546553, 0.43406308, 0.55474367))

# ------------------------------------------------------------------------------
# Export data for MicrobIEM
# ------------------------------------------------------------------------------

unique(otus_taxa$Genus)

taxonomy_microbIEM <- data.frame(
  Genus = c("Listeria", "Bacillus", "Pseudomonas", "Truepera", "Staphylococcus", 
            "Salmonella", "Enterococcus", "Imtechella", "Escherichia", 
            "Lactobacillus", "Allobacillus"),
  Taxonomy = c(
    "Bacteria;Firmicutes;Bacilli;Bacillales;Listeriaceae",
    "Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae",
    "Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae",
    "Bacteria;Deinococcus-Thermus;Deinococci;Deinococcales;Trueperaceae",
    "Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae",
    "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae",
    "Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae",
    "Bacteria;Bacteroidetes;Flavobacteria;Flavobacteriales;Flavobacteriaceae",
    "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae",
    "Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae",
    "Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae")
)
# Listeria: Bacteria;Firmicutes;Bacilli;Bacillales;Listeriaceae
# Bacillus: Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae
# Pseudomonas: Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae
# Truepera: Bacteria;Deinococcus-Thermus;Deinococci;Deinococcales;Trueperaceae
# Staphylococcus: Bacteria;Firmicutes;Bacilli;Bacillales;Staphylococcaceae
# Salmonella: Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae
# Enterococcus: Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae
# Imtechella: Bacteria;Bacteroidetes;Flavobacteria;Flavobacteriales;Flavobacteriaceae
# Escherichia: Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae
# Lactobacillus: Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae
# Allobacillus: Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae

otus_taxa_microbIEM <- merge(otus_taxa %>% rownames_to_column("OTU_ID"), 
                             taxonomy_microbIEM, by = "Genus", all = TRUE)


# Reformat otu table
otus_microbIEM <- otus %>% t() %>% as.data.frame() %>%
  merge(.,
        data.frame(
          OTU_ID = otus_taxa_microbIEM$OTU_ID,
          Taxonomy = with(otus_taxa_microbIEM, paste(Taxonomy, Genus, LV_4, 
                                                     sep = ";"))),
        by.x = 0, by.y = "OTU_ID") %>% 
  rename(OTU_ID = Row.names)

write.table(metadata, sep = "\t", row.names = F,
            paste0(file_directory, "Output/Data-formatted-for-MicrobIEM/MicrobIEM_mock_staggered-B_metafile.txt"))
write.table(otus_microbIEM, sep = "\t", row.names = F,
            paste0(file_directory, "Output/Data-formatted-for-MicrobIEM/MicrobIEM_mock_staggered-B_featurefile.txt"))
