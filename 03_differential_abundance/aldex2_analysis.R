################################################################################
# Differential Abundance Analysis with ALDEx2
# Description: Identify differentially abundant taxa between Sink A and Sink B
#              using ALDEx2 with IQLR normalization
################################################################################

# Load required libraries
library(ALDEx2)
library(dplyr)
library(ggplot2)

################################################################################
# DATA PREPARATION
################################################################################

# Load unrarefied data (ALDEx2 uses counts, not rarefied data)
feature_table_sinks <- read.table("physeq_feature_table_sinks_collapsed_unrar.txt", 
                                  header = TRUE, 
                                  row.names = 1, 
                                  sep = "\t")

metadata_sinks <- read.table("metadata_MSH.txt", 
                             header = TRUE, 
                             row.names = 1, 
                             sep = "\t")

# Remove counter samples
metadata_sinks <- metadata_sinks[metadata_sinks$Sample_Type != "Counter", ]

# Keep only Sink A and Sink B samples
metadata_sinks$Sink <- as.character(metadata_sinks$Sink)
metadata_sinks <- metadata_sinks[metadata_sinks$Sink %in% c("A", "B"), ]

# Remove 'X' prefix from sample names if present
colnames(feature_table_sinks) <- sub("^X", "", colnames(feature_table_sinks))

# Get sample names present in both datasets
common_samples <- intersect(colnames(feature_table_sinks), metadata_sinks$Sample_Name)

cat("=== Data Dimensions ===\n")
cat("Total samples in common:", length(common_samples), "\n")

################################################################################
# SUBSET AND ALIGN DATA
################################################################################

# Subset feature table to common samples
feature_table_filt <- feature_table_sinks[, common_samples]

# Subset and reorder metadata to match feature table columns
metadata_filt <- metadata_sinks[match(common_samples, metadata_sinks$Sample_Name), ]

# Verify alignment
if(!all(colnames(feature_table_filt) == metadata_filt$Sample_Name)) {
  stop("ERROR: Sample names don't match between feature table and metadata!")
}

cat("Feature table dimensions:", dim(feature_table_filt), "\n")
cat("Metadata rows:", nrow(metadata_filt), "\n")
cat("Samples match:", all(colnames(feature_table_filt) == metadata_filt$Sample_Name), "\n\n")

################################################################################
# FILTER LOW ABUNDANCE TAXA
################################################################################

# Calculate relative abundance
total_reads <- colSums(feature_table_filt)
rel_abundance <- sweep(feature_table_filt, 2, total_reads, FUN = "/")

# Remove taxa with maximum relative abundance < 1e-5 (0.001%)
# This reduces false positives
max_rel_abund <- apply(rel_abundance, 1, max)
features_to_keep <- max_rel_abund >= 1e-5

feature_table_filt <- feature_table_filt[features_to_keep, ]

cat("=== Filtering Low Abundance Taxa ===\n")
cat("Original number of features:", length(max_rel_abund), "\n")
cat("Features after filtering:", sum(features_to_keep), "\n")
cat("Features removed:", sum(!features_to_keep), "\n\n")

################################################################################
# RUN ALDEX2
################################################################################

# Create group vector
group <- metadata_filt$Sink

# Verify dimensions match
if(ncol(feature_table_filt) != length(group)) {
  stop("ERROR: Number of samples doesn't match group vector length!")
}

cat("=== Running ALDEx2 ===\n")
cat("Monte Carlo samples: 128\n")
cat("Test: Welch's t-test\n")
cat("Denominator: IQLR (inter-quartile log-ratio)\n\n")

# Run ALDEx2 analysis
# mc.samples = 128: Number of Monte Carlo instances
# test = "t": Welch's t-test
# effect = TRUE: Calculate effect sizes
# denom = "iqlr": Use inter-quartile log-ratio for normalization
set.seed(123)  # For reproducibility
aldex_results_sinks <- aldex(
  reads = feature_table_filt, 
  conditions = group, 
  mc.samples = 128, 
  test = "t", 
  effect = TRUE, 
  denom = "iqlr"
)

cat("ALDEx2 analysis complete!\n\n")

################################################################################
# ADD TAXONOMY INFORMATION
################################################################################

# Load taxonomy table
bacteria_taxonomy_unrar <- read.delim("physeq_taxonomy_sinks_collapsed_unrar.txt", 
                                      header = TRUE, 
                                      sep = '\t', 
                                      check.names = FALSE)
row.names(bacteria_taxonomy_unrar) <- bacteria_taxonomy_unrar$OTU_ID

# Merge ALDEx2 results with taxonomy
aldex_results_with_taxa <- merge(aldex_results_sinks, 
                                 bacteria_taxonomy_unrar, 
                                 by = "row.names", 
                                 all.x = TRUE)

# Rename first column
colnames(aldex_results_with_taxa)[1] <- "OTU_ID"

################################################################################
# IDENTIFY SIGNIFICANT TAXA
################################################################################

# Filter for significant taxa (p-value < 0.05)
significant_taxa <- aldex_results_with_taxa[aldex_results_with_taxa$we.ep < 0.05, ]

cat("=== Significant Taxa ===\n")
cat("Total significant taxa (p < 0.05):", nrow(significant_taxa), "\n")
cat("Enriched in Sink A (effect < 0):", sum(significant_taxa$effect < 0), "\n")
cat("Enriched in Sink B (effect > 0):", sum(significant_taxa$effect > 0), "\n\n")

# Get top taxa by effect size
top_taxa <- significant_taxa %>%
  dplyr::select(OTU_ID, Species, effect, we.ep, wi.ep, rab.all) %>%
  arrange(desc(abs(effect))) %>%
  head(50)

cat("=== Top 50 Differentially Abundant Taxa by Effect Size ===\n")
print(top_taxa)

################################################################################
# VISUALIZATION
################################################################################

# MA Plot (Mean abundance vs Difference)
ma_plot <- ggplot(aldex_results_sinks, 
                  aes(x = log2(rab.all + 1), y = diff.btw, 
                      color = we.ep < 0.05)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(
    values = c("TRUE" = "#d73027", "FALSE" = "grey70"),
    labels = c("TRUE" = "Significant (p < 0.05)", "FALSE" = "Not significant")
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(
    title = "ALDEx2 MA Plot",
    x = "Log2 Mean Abundance (rab.all)",
    y = "Difference Between Sinks (diff.btw)",
    color = "Significance"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("aldex2_ma_plot.pdf", ma_plot, width = 8, height = 6)

# Effect size plot
effect_plot <- ggplot(aldex_results_sinks, 
                      aes(x = log2(rab.all + 1), y = effect, 
                          color = we.ep < 0.05)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(
    values = c("TRUE" = "#d73027", "FALSE" = "grey70"),
    labels = c("TRUE" = "Significant (p < 0.05)", "FALSE" = "Not significant")
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted", color = "grey40") +
  theme_minimal() +
  labs(
    title = "ALDEx2 Effect Size Plot",
    x = "Log2 Mean Abundance (rab.all)",
    y = "Effect Size",
    color = "Significance"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("aldex2_effect_size_plot.pdf", effect_plot, width = 8, height = 6)

################################################################################
# SAVE RESULTS
################################################################################

# Save all ALDEx2 results
write.table(aldex_results_sinks, 
            "aldex2_results_sinks.tsv", 
            sep = "\t", 
            quote = FALSE, 
            col.names = NA)

# Save results with taxonomy
write.csv(aldex_results_with_taxa, 
          "aldex2_results_with_taxonomy.csv", 
          row.names = FALSE)

# Save significant taxa only
write.csv(significant_taxa, 
          "significant_da_taxa_sinks.csv", 
          row.names = FALSE)

# Save top taxa
write.csv(top_taxa, 
          "top_50_da_taxa_by_effect_size.csv", 
          row.names = FALSE)

cat("\n=== Files Saved ===\n")
cat("  - aldex2_results_sinks.tsv\n")
cat("  - aldex2_results_with_taxonomy.csv\n")
cat("  - significant_da_taxa_sinks.csv\n")
cat("  - top_50_da_taxa_by_effect_size.csv\n")
cat("  - aldex2_ma_plot.pdf\n")
cat("  - aldex2_effect_size_plot.pdf\n")

################################################################################
# NOTES:
# - ALDEx2 uses Monte Carlo sampling from Dirichlet distribution
# - IQLR (inter-quartile log-ratio) normalization reduces compositional bias
# - Effect size: difference between groups in CLR-transformed space
# - we.ep: Expected p-value from Welch's t-test
# - wi.ep: Expected p-value from Wilcoxon rank test
# - Negative effect: enriched in Sink A; Positive effect: enriched in Sink B
# - Taxa with abundance < 1e-5% removed to reduce false positives
################################################################################
