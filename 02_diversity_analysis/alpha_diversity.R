################################################################################
# Alpha Diversity Analysis
# Description: Chao1 and Shannon diversity analyses comparing sink communities
#              over time, using Wilcoxon and Kruskal-Wallis tests
################################################################################

# Load required libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(mctoolsr)

################################################################################
# DATA PREPARATION
################################################################################

# Import rarefied dataset (14,700 reads per sample)
bac <- load_taxa_table("asv_sinks_rar.txt", "sink_metadata_MSH_60days.txt")

# Remove control samples
bac <- filter_data(bac, 'Sample_Type_Specific', filter_vals = 'Control')
bac <- filter_data(bac, 'Sample_Type', 
                   filter_vals = c('Sterile P-trap Control', 
                                   'Filter Control', 
                                   'PCR Control'))

# Remove chloroplast, mitochondria, and unclassified taxa
bac <- filter_taxa_from_input(bac, 
                               taxa_to_remove = c('o__Chloroplast',
                                                  'f__Mitochondria',
                                                  'g__ unclassified'))

################################################################################
# CREATE PHYLOSEQ OBJECT
################################################################################

# Read metadata
bacteria_metadata_rar <- read.delim("sink_metadata_MSH_60days.txt", 
                                    header = TRUE, 
                                    sep = '\t', 
                                    row.names = 1)

# Read taxonomy
bacteria_taxonomy_rar <- read.delim(file = 'physeq_taxonomy_sinks_collapsed_rar.txt', 
                                    header = TRUE, 
                                    sep = '\t', 
                                    check.names = FALSE)
row.names(bacteria_taxonomy_rar) <- bacteria_taxonomy_rar$OTU_ID
bacteria_taxonomy_table_rar <- tax_table(as.matrix(bacteria_taxonomy_rar))

# Read ASV table
bacteria_asv_rar <- read.delim("physeq_feature_table_sinks_collapsed_rar.txt", 
                               header = TRUE,
                               sep = "\t", 
                               row.names = 1, 
                               check.names = FALSE)
bacteria_asv_table_rar <- otu_table(as.matrix(bacteria_asv_rar), taxa_are_rows = TRUE)

# Create phyloseq object
bacteria_physeq_rar <- phyloseq(
  otu_table(bacteria_asv_table_rar, taxa_are_rows = TRUE), 
  sample_data(bacteria_metadata_rar), 
  tax_table(bacteria_taxonomy_table_rar)
)

# Filter taxa and samples
bacteria_physeq_rar <- subset_taxa(bacteria_physeq_rar, 
                                   Order != "o__Chloroplast" &
                                   Family != "f__Mitochondria" &
                                   Genus != "g__ unclassified")

bacteria_physeq_rar <- subset_samples(bacteria_physeq_rar, Sink != "Control")

################################################################################
# CHAO1 DIVERSITY ANALYSIS
################################################################################

# Calculate Chao1 richness
chao <- estimate_richness(bacteria_physeq_rar, measures = c("Chao1"))
metadata_bac <- data.frame(sample_data(bacteria_physeq_rar))
chao_metadata <- cbind(metadata_bac, chao)

# Filter NA values and drop unused factor levels
chao_metadata <- chao_metadata %>%
  dplyr::filter(!is.na(Sink)) %>%
  mutate(Sink = droplevels(as.factor(Sink)))

chao_metadata$Week <- as.numeric(as.character(chao_metadata$Week))

# Overall comparison between sinks (Wilcoxon signed-rank test)
wilcox_result <- wilcox.test(Chao1 ~ Sink, data = chao_metadata)
cat("\n=== Overall Chao1 Comparison (Wilcoxon Test) ===\n")
print(wilcox_result)

# Summary statistics by sink and week
chao_summary <- chao_metadata %>%
  dplyr::group_by(Sink, Week) %>%
  dplyr::summarise(
    mean_chao = mean(Chao1, na.rm = TRUE),
    sd_chao = sd(Chao1, na.rm = TRUE),
    .groups = "drop"
  )

# Statistical tests by week (Wilcoxon for each week)
week_stats <- do.call(rbind, lapply(unique(chao_metadata$Week), function(w) {
  week_data <- chao_metadata[chao_metadata$Week == w, ]
  if(length(unique(week_data$Sink)) == 2) {
    p_val <- wilcox.test(Chao1 ~ Sink, data = week_data)$p.value
  } else {
    p_val <- NA_real_
  }
  data.frame(Week = as.numeric(w), p_value = p_val)
}))

# FDR correction (Bonferroni)
week_stats <- week_stats %>%
  mutate(
    p_adj = p.adjust(p_value, method = "fdr"),
    signif = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ NA_character_
    )
  )

cat("\n=== Chao1 by Week (Wilcoxon Tests with FDR Correction) ===\n")
print(week_stats)

# Calculate positions for significance stars
star_positions <- chao_summary %>%
  dplyr::group_by(Week) %>%
  dplyr::summarise(
    y_pos = max(mean_chao + sd_chao, na.rm = TRUE) * 1.1,
    .groups = "drop"
  ) %>%
  dplyr::mutate(Week = as.numeric(Week)) %>%
  dplyr::left_join(
    week_stats %>% dplyr::mutate(Week = as.numeric(Week)),
    by = "Week"
  ) %>%
  dplyr::filter(!is.na(signif))

# Create Chao1 plot
chao_plot_sd <- ggplot(chao_summary,
                       aes(x = Week, y = mean_chao,
                           color = Sink, group = Sink)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_chao - sd_chao,
                    ymax = mean_chao + sd_chao),
                width = 0.2) +
  geom_text(data = star_positions,
            aes(x = Week, y = y_pos, label = signif),
            color = "black",
            size = 5,
            inherit.aes = FALSE) +
  theme_minimal() +
  ylab("Mean Chao1 ± SD") +
  xlab("Week") +
  scale_x_continuous(breaks = sort(unique(chao_summary$Week))) +
  scale_color_manual(values = c("A" = "#4575b4", "B" = "#1a9850")) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  ) +
  labs(title = "Chao1 Diversity Over Time")

ggsave("chao1_diversity_by_week.pdf", chao_plot_sd, width = 8, height = 6)

################################################################################
# SHANNON DIVERSITY ANALYSIS
################################################################################

# Calculate Shannon diversity
shannon <- estimate_richness(bacteria_physeq_rar, measures = c("Shannon"))
shannon_metadata <- cbind(metadata_bac, shannon)

# Filter NA values
shannon_metadata <- shannon_metadata %>%
  dplyr::filter(!is.na(Sink)) %>%
  mutate(Sink = droplevels(as.factor(Sink)))

shannon_metadata$Week <- as.numeric(as.character(shannon_metadata$Week))

# Overall comparison between sinks (Wilcoxon test)
wilcox_shannon <- wilcox.test(Shannon ~ Sink, data = shannon_metadata)
cat("\n=== Overall Shannon Comparison (Wilcoxon Test) ===\n")
print(wilcox_shannon)

# Summary statistics
shannon_summary <- shannon_metadata %>%
  dplyr::group_by(Sink, Week) %>%
  dplyr::summarise(
    mean_shannon = mean(Shannon, na.rm = TRUE),
    sd_shannon = sd(Shannon, na.rm = TRUE),
    .groups = "drop"
  )

# Statistical tests by week
week_stats_shannon <- do.call(rbind, lapply(unique(shannon_metadata$Week), function(w) {
  week_data <- shannon_metadata[shannon_metadata$Week == w, ]
  if(length(unique(week_data$Sink)) == 2) {
    p_val <- wilcox.test(Shannon ~ Sink, data = week_data)$p.value
  } else {
    p_val <- NA_real_
  }
  data.frame(Week = as.numeric(w), p_value = p_val)
}))

week_stats_shannon <- week_stats_shannon %>%
  mutate(
    p_adj = p.adjust(p_value, method = "fdr"),
    signif = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ NA_character_
    )
  )

cat("\n=== Shannon by Week (Wilcoxon Tests with FDR Correction) ===\n")
print(week_stats_shannon)

# Create Shannon plot
star_positions_shannon <- shannon_summary %>%
  dplyr::group_by(Week) %>%
  dplyr::summarise(
    y_pos = max(mean_shannon + sd_shannon, na.rm = TRUE) * 1.05,
    .groups = "drop"
  ) %>%
  dplyr::mutate(Week = as.numeric(Week)) %>%
  dplyr::left_join(
    week_stats_shannon %>% dplyr::mutate(Week = as.numeric(Week)),
    by = "Week"
  ) %>%
  dplyr::filter(!is.na(signif))

shannon_plot_sd <- ggplot(shannon_summary,
                          aes(x = Week, y = mean_shannon,
                              color = Sink, group = Sink)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_shannon - sd_shannon,
                    ymax = mean_shannon + sd_shannon),
                width = 0.2) +
  geom_text(data = star_positions_shannon,
            aes(x = Week, y = y_pos, label = signif),
            color = "black",
            size = 5,
            inherit.aes = FALSE) +
  theme_minimal() +
  ylab("Mean Shannon ± SD") +
  xlab("Week") +
  scale_x_continuous(breaks = sort(unique(shannon_summary$Week))) +
  scale_color_manual(values = c("A" = "#4575b4", "B" = "#1a9850")) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  ) +
  labs(title = "Shannon Diversity Over Time")

ggsave("shannon_diversity_by_week.pdf", shannon_plot_sd, width = 8, height = 6)

################################################################################
# KRUSKAL-WALLIS AND DUNN'S POST-HOC TESTS
################################################################################

cat("\n========================================\n")
cat("KRUSKAL-WALLIS AND DUNN'S POST-HOC TESTS\n")
cat("========================================\n\n")

# Load required library for Dunn's test
library(dunn.test)

# Subset data by sink for within-sink temporal comparisons
sink_a_chao <- chao_metadata[chao_metadata$Sink == "A", ]
sink_b_chao <- chao_metadata[chao_metadata$Sink == "B", ]
sink_a_shannon <- shannon_metadata[shannon_metadata$Sink == "A", ]
sink_b_shannon <- shannon_metadata[shannon_metadata$Sink == "B", ]

# ============================================================================
# CHAO1 ANALYSIS
# ============================================================================

cat("--- Chao1 Diversity ---\n\n")

# Kruskal-Wallis tests for each sink across weeks
kw_chao_a <- kruskal.test(Chao1 ~ factor(Week), data = sink_a_chao)
kw_chao_b <- kruskal.test(Chao1 ~ factor(Week), data = sink_b_chao)

cat("Sink A - Kruskal-Wallis test:\n")
cat(sprintf("  Chi-squared = %.3f, df = %d, p-value = %.4f\n",
            kw_chao_a$statistic, kw_chao_a$parameter, kw_chao_a$p.value))

cat("\nSink B - Kruskal-Wallis test:\n")
cat(sprintf("  Chi-squared = %.3f, df = %d, p-value = %.4f\n\n",
            kw_chao_b$statistic, kw_chao_b$parameter, kw_chao_b$p.value))

# Dunn's post-hoc tests with Bonferroni correction
# (Only run when Kruskal-Wallis is significant, p < 0.05)

dunn_results_chao <- list()

if (kw_chao_a$p.value < 0.05) {
  cat("Sink A - Dunn's post-hoc test (Bonferroni correction):\n")
  
  dunn_chao_a <- dunn.test(sink_a_chao$Chao1,
                           factor(sink_a_chao$Week),
                           method = "bonferroni",
                           kw = FALSE,
                           table = FALSE)
  
  # Extract significant pairwise comparisons
  sig_idx <- which(dunn_chao_a$P.adjusted < 0.05)
  
  if (length(sig_idx) > 0) {
    sig_comparisons_a <- data.frame(
      Sink = "A",
      Metric = "Chao1",
      Comparison = dunn_chao_a$comparisons[sig_idx],
      Z = dunn_chao_a$Z[sig_idx],
      P.adjusted = dunn_chao_a$P.adjusted[sig_idx]
    )
    print(sig_comparisons_a)
    dunn_results_chao$sink_a <- sig_comparisons_a
  } else {
    cat("  No significant pairwise differences after Bonferroni correction\n")
  }
  cat("\n")
} else {
  cat("Sink A - Kruskal-Wallis not significant; Dunn's test not performed\n\n")
}

if (kw_chao_b$p.value < 0.05) {
  cat("Sink B - Dunn's post-hoc test (Bonferroni correction):\n")
  
  dunn_chao_b <- dunn.test(sink_b_chao$Chao1,
                           factor(sink_b_chao$Week),
                           method = "bonferroni",
                           kw = FALSE,
                           table = FALSE)
  
  # Extract significant pairwise comparisons
  sig_idx <- which(dunn_chao_b$P.adjusted < 0.05)
  
  if (length(sig_idx) > 0) {
    sig_comparisons_b <- data.frame(
      Sink = "B",
      Metric = "Chao1",
      Comparison = dunn_chao_b$comparisons[sig_idx],
      Z = dunn_chao_b$Z[sig_idx],
      P.adjusted = dunn_chao_b$P.adjusted[sig_idx]
    )
    print(sig_comparisons_b)
    dunn_results_chao$sink_b <- sig_comparisons_b
  } else {
    cat("  No significant pairwise differences after Bonferroni correction\n")
  }
  cat("\n")
} else {
  cat("Sink B - Kruskal-Wallis not significant; Dunn's test not performed\n\n")
}

# ============================================================================
# SHANNON ANALYSIS
# ============================================================================

cat("--- Shannon Diversity ---\n\n")

# Kruskal-Wallis tests for each sink across weeks
kw_shannon_a <- kruskal.test(Shannon ~ factor(Week), data = sink_a_shannon)
kw_shannon_b <- kruskal.test(Shannon ~ factor(Week), data = sink_b_shannon)

cat("Sink A - Kruskal-Wallis test:\n")
cat(sprintf("  Chi-squared = %.3f, df = %d, p-value = %.4f\n",
            kw_shannon_a$statistic, kw_shannon_a$parameter, kw_shannon_a$p.value))

cat("\nSink B - Kruskal-Wallis test:\n")
cat(sprintf("  Chi-squared = %.3f, df = %d, p-value = %.4f\n\n",
            kw_shannon_b$statistic, kw_shannon_b$parameter, kw_shannon_b$p.value))

dunn_results_shannon <- list()

if (kw_shannon_a$p.value < 0.05) {
  cat("Sink A - Dunn's post-hoc test (Bonferroni correction):\n")
  
  dunn_shannon_a <- dunn.test(sink_a_shannon$Shannon,
                              factor(sink_a_shannon$Week),
                              method = "bonferroni",
                              kw = FALSE,
                              table = FALSE)
  
  # Extract significant pairwise comparisons
  sig_idx <- which(dunn_shannon_a$P.adjusted < 0.05)
  
  if (length(sig_idx) > 0) {
    sig_comparisons_a <- data.frame(
      Sink = "A",
      Metric = "Shannon",
      Comparison = dunn_shannon_a$comparisons[sig_idx],
      Z = dunn_shannon_a$Z[sig_idx],
      P.adjusted = dunn_shannon_a$P.adjusted[sig_idx]
    )
    print(sig_comparisons_a)
    dunn_results_shannon$sink_a <- sig_comparisons_a
  } else {
    cat("  No significant pairwise differences after Bonferroni correction\n")
  }
  cat("\n")
} else {
  cat("Sink A - Kruskal-Wallis not significant; Dunn's test not performed\n\n")
}

if (kw_shannon_b$p.value < 0.05) {
  cat("Sink B - Dunn's post-hoc test (Bonferroni correction):\n")
  
  dunn_shannon_b <- dunn.test(sink_b_shannon$Shannon,
                              factor(sink_b_shannon$Week),
                              method = "bonferroni",
                              kw = FALSE,
                              table = FALSE)
  
  # Extract significant pairwise comparisons
  sig_idx <- which(dunn_shannon_b$P.adjusted < 0.05)
  
  if (length(sig_idx) > 0) {
    sig_comparisons_b <- data.frame(
      Sink = "B",
      Metric = "Shannon",
      Comparison = dunn_shannon_b$comparisons[sig_idx],
      Z = dunn_shannon_b$Z[sig_idx],
      P.adjusted = dunn_shannon_b$P.adjusted[sig_idx]
    )
    print(sig_comparisons_b)
    dunn_results_shannon$sink_b <- sig_comparisons_b
  } else {
    cat("  No significant pairwise differences after Bonferroni correction\n")
  }
  cat("\n")
} else {
  cat("Sink B - Kruskal-Wallis not significant; Dunn's test not performed\n\n")
}

# Combine all significant Dunn's test results
all_dunn_results <- do.call(rbind, c(dunn_results_chao, dunn_results_shannon))

################################################################################
# SAVE RESULTS
################################################################################

# Save results
write.csv(chao_summary, "chao1_summary_statistics.csv", row.names = FALSE)
write.csv(shannon_summary, "shannon_summary_statistics.csv", row.names = FALSE)
write.csv(week_stats, "chao1_weekly_tests.csv", row.names = FALSE)
write.csv(week_stats_shannon, "shannon_weekly_tests.csv", row.names = FALSE)

# Save Kruskal-Wallis results
kw_results <- data.frame(
  Sink = c("A", "B", "A", "B"),
  Metric = c("Chao1", "Chao1", "Shannon", "Shannon"),
  Chi_squared = c(kw_chao_a$statistic, kw_chao_b$statistic,
                  kw_shannon_a$statistic, kw_shannon_b$statistic),
  df = c(kw_chao_a$parameter, kw_chao_b$parameter,
         kw_shannon_a$parameter, kw_shannon_b$parameter),
  p_value = c(kw_chao_a$p.value, kw_chao_b$p.value,
              kw_shannon_a$p.value, kw_shannon_b$p.value)
)
write.csv(kw_results, "kruskal_wallis_results.csv", row.names = FALSE)

# Save Dunn's test results if any significant comparisons exist
if (!is.null(all_dunn_results) && nrow(all_dunn_results) > 0) {
  write.csv(all_dunn_results, "dunn_posthoc_significant_results.csv", row.names = FALSE)
  cat("\nSignificant Dunn's post-hoc results saved to dunn_posthoc_significant_results.csv\n")
}

cat("\n=== Alpha Diversity Analysis Complete ===\n")
cat("\nFiles saved:\n")
cat("  - chao1_summary_statistics.csv\n")
cat("  - shannon_summary_statistics.csv\n")
cat("  - chao1_weekly_tests.csv\n")
cat("  - shannon_weekly_tests.csv\n")
cat("  - kruskal_wallis_results.csv\n")
cat("  - dunn_posthoc_significant_results.csv (if significant results found)\n")

################################################################################
# NOTES:
# - Wilcoxon signed-rank test used for overall sink comparisons (paired data)
# - Weekly comparisons use Wilcoxon tests with FDR correction for visualization
# - Kruskal-Wallis test evaluates temporal variation within each sink
# - Dunn's post-hoc test (Bonferroni-corrected) identifies specific week-to-week
#   differences when Kruskal-Wallis is significant
# - Rarefied to 14,700 reads per sample (n=9 samples lost)
################################################################################
