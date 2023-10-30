
# ------------------------------------------------------------------------------
# Run frequency filter
# ------------------------------------------------------------------------------

source(paste0(file_directory, "BM_ST_PR.R"))
res_frequency_staggered_PR <- BM_frequency_PR(
  otus_rel = otus_rel, metadata = metadata, Mock_info = Tax_class, 
  dmax = 7) # spike: 1, dmock: 7

# ------------------------------------------------------------------------------
# Run MicrobIEM ratio filter
# ------------------------------------------------------------------------------

res_microbiem_ratio_staggered_PR <- BM_microbiem_ratio_PR(
  otus_rel = otus_rel, metadata = metadata, Mock_info = Tax_class, 
  dmax = 7) # spike: 1, dmock: 7

# ------------------------------------------------------------------------------
# Run MicrobIEM span filter
# ------------------------------------------------------------------------------

res_microbiem_span_staggered_PR <- BM_microbiem_span_PR(
  otus_rel = otus_rel, metadata = metadata, Mock_info = Tax_class, 
  dmax = 7) # spike: 1, dmock: 7

# ------------------------------------------------------------------------------
# Run Decontam prevalence filter
# ------------------------------------------------------------------------------

res_decontprev_staggered_PR <- BM_decontam_prev_PR(
  otus_rel = otus_rel, metadata = metadata, Mock_info = Tax_class, 
  per_dilution = TRUE, dmax = 7) # spike: 1, dmock: 7

# ------------------------------------------------------------------------------
# Run Decontam frequency filter
# ------------------------------------------------------------------------------

res_decontfreq_staggered_PR <- BM_decontam_freq_PR(
  otus_rel = otus_rel, metadata = metadata, Mock_info = Tax_class, 
  per_dilution = TRUE, dmax = 7) # spike: 1, dmock: 7

# ------------------------------------------------------------------------------
# Run SourceTracker
# ------------------------------------------------------------------------------

#source(paste0(file_directory, "BM_ST_PR.R"))
print(Sys.time()) # takes 2.5-3 hours
res_SourceTracker_staggered_PR <- BM_sourcetracker_PR(
  file_directory = file_directory,
  otus = otus, metadata = metadata, Mock_info = Tax_class, 
  raref_depth = 1000) # spike: 1, dmock: 7
print(Sys.time())
rm(BM_sourcetracker_PR)

# Save all:
###save(list = ls()[grepl("res_.*_PR", ls())], file = paste0(file_directory, "Output/R_objects/2_BM_staggered_all_PR.RData"))
###load(paste0(file_directory, "Output/R_objects/2_BM_staggered_all_PR.RData"), verbose = TRUE)

################################################################################

# Combine all benchmarking results
combined_results_PR <- Reduce(
  rbind.data.frame, list(
    data.frame(res_frequency_staggered_PR, Method = "Frequency"),
    data.frame(res_microbiem_ratio_staggered_PR, Method = "MicrobIEMRatio"),
    data.frame(res_microbiem_span_staggered_PR, Method = "MicrobIEMSpan") %>%
      mutate(Filter = substring(Filter, 1, 22)),
    data.frame(res_decontfreq_staggered_PR, Method = "DecontamFreq"), 
    data.frame(filter(res_SourceTracker_staggered_PR, grepl(", a1", res_SourceTracker_staggered_PR$Filter)), Method = "SourceTracker, a1"),
    data.frame(filter(res_SourceTracker_staggered_PR, grepl(", a2", res_SourceTracker_staggered_PR$Filter)), Method = "SourceTracker, a2"),
    data.frame(filter(res_SourceTracker_staggered_PR, grepl(", b1", res_SourceTracker_staggered_PR$Filter)), Method = "SourceTracker, b1"),
    data.frame(res_decontprev_staggered_PR, Method = "DecontamPrev")))

# Calculate Accuracy, Sensitivity, Specificity, Youden's index, Matthews index
combined_results_PR <- combined_results_PR %>%
  mutate(Accuracy = (TP + TN) / (TP + TN + FP + FN),
         Sensitivity = TP / (TP+FN),
         Precision = TP / (TP+FP),
         Specificity = TN / (TN+FP),
         Youden = Sensitivity + Specificity - 1,
         Matthews = (TP*TN - FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))

#combined_results_PR <- combined_results_PR %>% 
#  mutate(Filter = gsub(", a2", ", a2 = ", Filter),
#         Filter = gsub(", a1", ", a1 = ", Filter),
#         Filter = gsub(", b1", ", b1 = ", Filter))
combined_results_PR <- combined_results_PR %>%
  mutate(Precision = case_when(
    (is.na(Precision) & TP == 0 & FP == 0 & TN > 0 & FN > 0) ~ 0,
    TRUE ~ Precision))
combined_results_PR <- combined_results_PR %>% 
  mutate(Method = factor(Method, levels = c(
  "Frequency", sort(unique(
    combined_results_PR$Method[combined_results_PR$Method != "Frequency"])))))

combined_results_PR %>% 
  merge(., metadata, by.x = "Sample", by.y = "Sample_ID") %>%
  filter(Dilution %in% c("D1", "D2", "D3", "D4")) %>% # only in Dmock
  group_by(Dilution, Filter, Method) %>%
  summarise(Sensitivity = mean(Sensitivity, na.rm = T),
            Precision = mean(Precision, na.rm= T)) %>% #View
  mutate(Filter = gsub(".* = ", "", Filter, perl = T)) %>% 
  ggplot(., aes(x = Sensitivity, y = Precision, label = Filter, colour = Dilution)) +
  geom_point() + geom_line(aes(group = Dilution)) + 
  geom_text(hjust = -0.2, vjust = 0, size = 2) +
  facet_wrap(. ~ Method, nrow = 2,
             labeller = labeller(Method = c("Frequency" = "Frequency filter",
                                            "DecontamFreq" = "Decontam (freq.)",
                                            "DecontamPrev" = "Decontam (prev.)",
                                            "PresenceNEG2" = "Presence filter",
                                            "MicrobIEMRatio" = "MicrobIEM (ratio)",
                                            "MicrobIEMSpan" = "MicrobIEM (span)",
                                            "SourceTracker, a1" = "SourceTracker (a1)",
                                            "SourceTracker, a2" = "SourceTracker (a2)",
                                            "SourceTracker, b1" = "SourceTracker (b1)"))) +
  ggtitle("Staggered mock community A") +
  scale_colour_manual(values = theme_colours[c(1, 3, 11, 5)]) +
  plot_theme + theme(plot.title = element_text(size = rel(1.1))) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"))
ggsave(paste0(file_directory, "Output/Plots/PR-curves_dmock.svg"),
       width = 8.2, height = 4.5)

