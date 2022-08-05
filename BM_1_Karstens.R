################################################################################
#
# Decontamination benchmarking on an even mock community (Karstens data) 
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

# Download Karstens data from Github:
# https://github.com/lakarstens/ControllingContaminants16S/blob/master/Analyses/mockDilutionsPrep.RData
# Load the R object that contains all data
load(paste0(file_directory, "Input/Karstens/mockDilutionsPrep.RData"), 
     verbose = TRUE)

# Define otu table with count data (samples in rows, ASVs in columns)
otus <- as.data.frame(ps@otu_table@.Data)

# Create OTU table with relative abundances
otus_rel <- sweep(otus, 1, rowSums(otus), "/")
all(rowSums(otus_rel) == 1) # Check

# Create metadata
metadata <- data.frame(
   Sample_ID = row.names(otus),
   Sample_type = ifelse(row.names(otus) == "Blank", "NEG2", "SAMPLE")
)

# Add DNA concentration for decontam
metadata <- merge(metadata, mock_ps@sam_data[, "DNA_conc"], by.x = "Sample_ID",
                  by.y = 0, all = TRUE)
# Define number of input cells
metadata["Input_cells"] <- c(0, 1.5e+9/3^c(0:8))
# Define the dilution
metadata["Dilution"] <- ifelse(metadata$Sample_ID == "Blank", "NEG", 
                               metadata$Sample_ID)
# Add Source/sink information for SourceTracker
metadata["SourceSink"] <- ifelse(metadata$Sample_type == "NEG2", "Source", 
                                 "Sink")

# Create gold standard (classify ASVs into mock/contaminant by sequence)
Tax_class <- data.frame(
   Contaminant = ifelse(colnames(otus) %in% mock_taxa, "Mock", "Contaminant"),
   row.names = colnames(otus))

# Save taxonomic information
otus_taxa <- as.data.frame(ps@tax_table@.Data)
# Karstens et al. compared the ASVs sequences with the expected sequences, 
#  and show "ASV_8" as "Salmonella_ASV2" in their plots. 
# Therefore, we correct the taxonomic assignment of ASV_8 to Salmonella.
otus_taxa["ASV_8", "Genus"] <- "Salmonella"
# "Lactobacillus fermentum" is now called "Limosilactobacillus fermentum".
otus_taxa["ASV_1", "Genus"] <- "Limosilactobacillus"

# Remove all objects from the workspace that are not needed anymore
rm(list = ls()[!ls() %in% c("file_directory", "otus", "otus_rel", "otus_taxa", 
                            "metadata", "Tax_class", "plot_theme", 
                            "theme_colours", "benchm_colours")])

# ------------------------------------------------------------------------------
# Plot sample composition per dilution
# ------------------------------------------------------------------------------

# https://stackoverflow.com/questions/63010394/superscript-in-axis-labels-in-ggplot2-for-ions
svg("Output/Plots/Taxonomy_even.svg", width = 6, height = 3.8)
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
   merge(., metadata, by.x = "name", by.y = "Sample_ID", all.x = TRUE) %>%
   # Add fake row to separate samples and the control by an empty bar
   add_row(Category = "Contaminants", Dilution = "", value = NA) %>%
   # Control order of genera and samples
   mutate(Category = factor(Category, levels = c(unique(sort(.$Category))[c(1, 3:9)], "Contaminants")),
      Dilution = factor(Dilution, levels = c(paste0("D", 0:8), "", "NEG"))) %>%
   ggplot(aes(x = Dilution, y = value, fill = Category)) +
   geom_bar(stat = "identity") +
   #xlab(bquote('1:3 dilution series (1.5 x'~10^9*' to 2.3 x'~10^5*')')) + 
   xlab(bquote('1:3 dilution series ('*10^9*' to'~10^5*' cells)')) + 
   ylab("Relative abundance") +
   plot_theme +
   scale_y_continuous(expand = c(0, 0)) +
   theme(axis.ticks.x = element_blank(),
         axis.text.x = element_text(face = c("plain", rep(c("plain", "bold"), 4), rep("plain", 2)),
                                    colour = c("grey30", rep(c("grey30", "black"), 4), rep("grey30", 2)))) +
   scale_fill_manual(
      "Genus", values = c(theme_colours[c(20, 21, 9, 18, 14, 5, 7, 23)], "grey80"))
dev.off()

# ------------------------------------------------------------------------------
# Run SourceTracker
# ------------------------------------------------------------------------------

source(paste0(file_directory, "BM_SourceTracker.R"))
res_sourcetracker_karstens <- BM_sourcetracker(
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
res_frequency_karstens <- BM_frequency(
   otus_rel = otus_rel, metadata = metadata, Mock_info = Tax_class, dmax = 8)
rm(BM_frequency)

# ------------------------------------------------------------------------------
# Run NEG2 presence filter
# ------------------------------------------------------------------------------

source(paste0(file_directory, "BM_NegPresence.R"))
res_negpresence_karstens <- BM_negpresence(
   otus_rel = otus_rel, metadata = metadata, Mock_info = Tax_class, dmax = 8)
rm(BM_negpresence)

# ------------------------------------------------------------------------------
# Run MicrobIEM ratio filter
# ------------------------------------------------------------------------------

source(paste0(file_directory, "BM_MicrobIEM.R"))
res_microbiem_karstens <- BM_microbiem(
   otus_rel = otus_rel, metadata = metadata, Mock_info = Tax_class, dmax = 8)
rm(MicrobIEM_decontamination, BM_microbiem)

# ------------------------------------------------------------------------------
# Run Decontam prevalence filter
# ------------------------------------------------------------------------------

source(paste0(file_directory, "BM_DecontamPrev.R"))
res_decontprev_karstens <- BM_decontam_prev(
   otus_rel = otus_rel, metadata = metadata, Mock_info = Tax_class, 
   per_dilution = TRUE, dmax = 8)
rm(BM_decontam_prev)

# ------------------------------------------------------------------------------
# Run Decontam frequency filter
# ------------------------------------------------------------------------------

source(paste0(file_directory, "BM_DecontamFreq.R"))
res_decontfreq_karstens <- BM_decontam_freq(
   otus_rel = otus_rel, metadata = metadata, Mock_info = Tax_class, 
   per_dilution = FALSE, dmax = 8)
rm(BM_decontam_freq)

# ------------------------------------------------------------------------------
# Save the results
# ------------------------------------------------------------------------------

###save(list = ls()[grepl("res_", ls())], file = paste0(file_directory, "Output/R_objects/Karstens_benchm_res.RData"))
###load(paste0(file_directory, "Output/R_objects/Karstens_benchm_res.RData"), verbose = TRUE)

################################################################################
#
# Combination of results 
#
################################################################################

# Combine all benchmarking results
combined_results <- Reduce(
   rbind.data.frame, list(
      data.frame(res_sourcetracker_karstens, Method = "SourceTracker"),
      data.frame(res_frequency_karstens, Method = "Frequency"),
      data.frame(res_negpresence_karstens, Method = "PresenceNEG2"),
      data.frame(res_microbiem_karstens, Method = "MicrobIEM"),
      data.frame(res_decontfreq_karstens, Method = "DecontamFreq"), 
      data.frame(res_decontprev_karstens, Method = "DecontamPrev")))

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
combined_results["Dilution"] <- combined_results$Sample

# Order methods, filters, and measures for plots
combined_results$Method <- 
   factor(combined_results$Method, levels = c(
      "Frequency", "DecontamFreq", "DecontamPrev", "PresenceNEG2", 
      "SourceTracker", "MicrobIEM"))
combined_results$Measure <- 
   factor(combined_results$Measure, levels = c(
      "Accuracy", "Sensitivity", "Specificity", "Youden", "Matthews"))
combined_results$Filter <- factor(combined_results$Filter, levels = c(
   unique(res_frequency_karstens$Filter),
   unique(res_decontfreq_karstens$Filter),
   unique(res_decontprev_karstens$Filter),
   unique(res_negpresence_karstens$Filter),
   unique(res_sourcetracker_karstens$Filter),
   unique(res_microbiem_karstens$Filter)))

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
   t() %>% as.data.frame()
names(Cont_prev) <- "Cont_prev"
# Add contaminant prevalence to combined results
combined_results <- merge(combined_results, Cont_prev, by.x = "Sample", 
                          by.y = 0, all.x = TRUE)
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

svg(paste0(file_directory, "Output/Plots/Benchm_even_supplement.svg"), 
    width = 9.1, height = 6.5) # width = 12, height = 6
p <- 
   combined_results %>%
   # Calculate mean per dilution (not needed, only for compatibility with Dmock)
   group_by(Filter, Dilution, Measure) %>%
   mutate(Mean = mean(Value)) %>% 
   # Reduce shown filters & select samples
   filter(!Filter %in% c("MicrobIEM, ratio = 1.5", "SourceTracker, beta = 0.001",
                         "Decontam (frequency), thr = 1", 
                         "Decontam (prevalence), thr = 1"),
          Dilution %in% c("D2", "D4", "D6", "D8")) %>% 
   ggplot(aes(x = Cont_prev, y = Mean, colour = Filter, group = Filter)) +
   plot_theme +
   theme(panel.grid.major = element_line(colour = "#e0e0e0"),
         legend.position = "none") +
   # Add base lines for Accuracy and Youden's index
   geom_abline(aes(slope = -Cont_prev_line, intercept = 1), colour = "black") +
   geom_abline(aes(slope = Youden_line, intercept = 0), colour = "black") +
   xlab("Contaminant prevalence") +
   ylab("Value per dilution") +
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

svg(paste0(file_directory, "Output/Plots/Benchm_even_main.svg"), 
    width = 8, height = 3.5) # height = 3.31
p <- 
   combined_results %>%
   group_by(Method, Filter, Dilution, Measure, X_axis_prox) %>%
   summarise(Mean = mean(Value)) %>%
   filter(Measure %in% c("Youden"),
          !Filter %in% c("MicrobIEM, ratio = 1.5", "SourceTracker, beta = 0.001",
                         "Decontam (frequency), thr = 1", 
                         "Decontam (prevalence), thr = 1"),
          Dilution %in% c("D2", "D4", "D6", "D8")) %>% 
   as.data.frame() %>%
   mutate(Dilution = factor(Dilution, labels = c("D2~(1.7~x~10^8)",
                                                 "D4~(1.9~x~10^7)",
                                                 "D6~(2.1~x~10^6)",
                                                 "D8~(2.3~x~10^5)"))) %>% 
   ggplot(aes(x = X_axis_prox, y = Mean, colour = Filter, fill = Filter)) +
   plot_theme +
   theme(panel.grid.major.y = element_line(colour = "#e0e0e0"),
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
# Legend
# ------------------------------------------------------------------------------

for (i in c("MicrobIEM", "SourceTracker", "PresenceNEG2", "DecontamPrev",
            "DecontamFreq", "Frequency"
            )) {
   p <- 
      combined_results %>%
      filter(Method == i) %>% 
      filter(!Filter %in% c("MicrobIEM, ratio = 1.5", "SourceTracker, beta = 0.001",
                            "Decontam (frequency), thr = 1", 
                            "Decontam (prevalence), thr = 1"),
             Dilution %in% c("D2", "D4", "D6", "D8")
             ) %>% as.data.frame() %>%
      mutate(Filter = gsub(".*, ", "", Filter)) %>% 
      {if (i %in% c("MicrobIEM", "SourceTracker")) 
         mutate(., Filter = factor(Filter, levels = gsub(
            ".*, ", "", levels(combined_results$Filter))[c(15:17, 19, 21:23)])) else .} %>%
      ggplot(aes(x = X_axis_prox, y = Value, colour = Filter, fill = Filter)) +
      plot_theme +
      theme(legend.spacing.y = unit(-0.1, 'cm'),
            legend.position = "left",
            legend.title = element_text(margin = margin(0, 0, 6, 0)),
            legend.justification = c(0, 0.8)) +
      geom_col(position = position_dodge2(preserve = "single")) +
      guides(colour = guide_legend(ncol = 1, byrow = TRUE))
      if (i == "MicrobIEM") {
         p1 <- p + scale_colour_manual(name = "MicrobIEM", values = benchm_colours[16:19]) +
            scale_fill_manual(name = "MicrobIEM", values = benchm_colours[16:19])
      }
      if (i == "SourceTracker") {
         p2 <- p + scale_colour_manual(name = "SourceTracker", values = benchm_colours[13:15]) +
            scale_fill_manual(name = "SourceTracker", values = benchm_colours[13:15])
      }
      if (i == "PresenceNEG2") {
         p3 <- p + scale_colour_manual(name = "Presence filter", values = benchm_colours[12]) +
            scale_fill_manual(name = "Presence filter", values = benchm_colours[12])
      }
      if (i == "DecontamPrev") {
         p4 <- p + scale_colour_manual(name = "Decontam (prevalence)", values = benchm_colours[8:11]) +
            scale_fill_manual(name = "Decontam (prevalence)", values = benchm_colours[8:11])
      }
      if (i == "DecontamFreq") {
         p5 <- p + scale_colour_manual(name = "Decontam (frequency)", values = benchm_colours[4:7]) +
            scale_fill_manual(name = "Decontam (frequency)", values = benchm_colours[4:7])
      }
      if (i == "Frequency") {
         p6 <- p + scale_colour_manual(name = "Frequency filter", values = benchm_colours) +
            scale_fill_manual(name = "Frequency filter", values = benchm_colours)
      }
}

svg(paste0(file_directory, "Output/Plots/Benchm_legend.svg"),
    width = 1.8, height = 6.1)
plot_grid(as_ggplot(get_legend(p6)), as_ggplot(get_legend(p5)),
          as_ggplot(get_legend(p4)), as_ggplot(get_legend(p3)), 
          as_ggplot(get_legend(p2)), as_ggplot(get_legend(p1)), ncol = 1, 
          align = "hv", rel_heights = c(0.91, 1.1, 1.1, 0.55, 0.91, 1.1))
dev.off()

# ------------------------------------------------------------------------------
# Text data
# ------------------------------------------------------------------------------

# Contaminant prevalence per sample and NEG
merge(t(otus_rel), Tax_class, by = 0) %>% 
   group_by(Contaminant) %>%
   summarise_at(.vars = c(Sample_IDs, "Blank"), .funs = sum) %>% 
   filter(Contaminant == "Contaminant") %>% 
   select(all_of(c(Sample_IDs, "Blank"))) %>% 
   t() %>% as.data.frame()
#D0    0.0005011554
#D1    0.0013879652
#D2    0.0179780835
#D3    0.0445470371
#D4    0.1195566259
#D5    0.2786355186
#D6    0.6448189422
#D7    0.5576927759
#D8    0.8006450509
#Blank 0.9882705673
