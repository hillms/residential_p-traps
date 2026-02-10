################################################################################
# Beta Diversity Analysis
# Description: Community composition analysis using Aitchison distance and
#              PERMANOVA to compare sink bacterial communities
################################################################################

# Load required libraries
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(compositions)

################################################################################
# DATA PREPARATION
################################################################################

# Load phyloseq object (created in alpha_diversity.R)
# bacteria_physeq_rar should already be loaded

# Extract components
metadata <- data.frame(sample_data(bacteria_physeq_rar))
otu_table_mat <- as.matrix(otu_table(bacteria_physeq_rar))

################################################################################
# CENTERED LOG-RATIO (CLR) TRANSFORMATION
################################################################################

# Transpose OTU table (samples as rows, ASVs as columns)
otu_t <- t(otu_table_mat)

# Add pseudocount to handle zeros (required for CLR transformation)
otu_t_pseudo <- otu_t + 0.5

# Apply CLR transformation using compositions package
otu_clr <- clr(otu_t_pseudo)

cat("\n=== CLR Transformation Complete ===\n")
cat("Transformed table dimensions:", dim(otu_clr), "\n")

################################################################################
# AITCHISON DISTANCE CALCULATION
################################################################################

# Calculate Aitchison distance (Euclidean distance on CLR-transformed data)
aitchison_dist <- dist(otu_clr, method = "euclidean")

cat("\n=== Aitchison Distance Matrix Calculated ===\n")
cat("Distance matrix size:", attr(aitchison_dist, "Size"), "samples\n")

################################################################################
# PERMANOVA ANALYSIS
################################################################################

# Overall PERMANOVA - comparing Sink A vs Sink B
set.seed(123)  # For reproducibility
permanova_sink <- adonis2(aitchison_dist ~ Sink, 
                         data = metadata, 
                         permutations = 999,
                         method = "euclidean")

cat("\n=== PERMANOVA: Sink A vs Sink B ===\n")
print(permanova_sink)

# PERMANOVA by Week
# Test if communities differ between sinks within each week
week_permanova_results <- data.frame()

for(week in unique(metadata$Week)) {
  # Subset data for this week
  week_samples <- rownames(metadata[metadata$Week == week, ])
  week_indices <- which(rownames(metadata) %in% week_samples)
  
  if(length(week_indices) > 3) {  # Need sufficient samples
    week_dist <- as.dist(as.matrix(aitchison_dist)[week_indices, week_indices])
    week_meta <- metadata[week_indices, ]
    
    # Run PERMANOVA for this week
    set.seed(123)
    perm_result <- adonis2(week_dist ~ Sink, 
                          data = week_meta, 
                          permutations = 999,
                          method = "euclidean")
    
    week_permanova_results <- rbind(week_permanova_results,
                                   data.frame(
                                     Week = week,
                                     R2 = perm_result$R2[1],
                                     F_statistic = perm_result$F[1],
                                     p_value = perm_result$`Pr(>F)`[1]
                                   ))
  }
}

# FDR correction for multiple testing
week_permanova_results$p_adj <- p.adjust(week_permanova_results$p_value, method = "fdr")

cat("\n=== PERMANOVA Results by Week ===\n")
print(week_permanova_results)

################################################################################
# PRINCIPAL COORDINATES ANALYSIS (PCoA)
################################################################################

# Perform PCoA on Aitchison distances
pcoa_result <- cmdscale(aitchison_dist, k = 2, eig = TRUE)

# Extract variance explained
variance_explained <- pcoa_result$eig / sum(pcoa_result$eig) * 100

# Create PCoA data frame
pcoa_df <- data.frame(
  PC1 = pcoa_result$points[, 1],
  PC2 = pcoa_result$points[, 2],
  Sink = metadata$Sink,
  Week = metadata$Week
)

cat("\n=== PCoA Variance Explained ===\n")
cat("PC1:", round(variance_explained[1], 2), "%\n")
cat("PC2:", round(variance_explained[2], 2), "%\n")

# Create PCoA plot
pcoa_plot <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Sink, shape = Sink)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  scale_color_manual(values = c("A" = "#4575b4", "B" = "#1a9850")) +
  labs(
    x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
    title = "PCoA of Sink Communities (Aitchison Distance)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("pcoa_aitchison_distance.pdf", pcoa_plot, width = 8, height = 6)

# PCoA colored by week
pcoa_plot_week <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = as.factor(Week))) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_viridis_d() +
  labs(
    x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
    title = "PCoA of Sink Communities by Week",
    color = "Week"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave("pcoa_aitchison_by_week.pdf", pcoa_plot_week, width = 10, height = 6)

################################################################################
# BRAY-CURTIS DISSIMILARITY (for temporal stability analysis)
################################################################################

# Calculate Bray-Curtis dissimilarity on rarefied counts
bray_dist <- vegdist(otu_t, method = "bray")

# PERMANOVA with Bray-Curtis
set.seed(123)
permanova_bray <- adonis2(bray_dist ~ Sink, 
                         data = metadata, 
                         permutations = 999)

cat("\n=== PERMANOVA: Bray-Curtis Dissimilarity ===\n")
print(permanova_bray)

# Calculate beta diversity from baseline (timepoint 1) for each sink
# This is used in temporal stability analysis

beta_div_from_baseline <- function(dist_matrix, metadata, baseline_timepoint = 1) {
  results <- data.frame()
  
  for(sink_id in unique(metadata$Sink)) {
    sink_meta <- metadata[metadata$Sink == sink_id, ]
    sink_samples <- rownames(sink_meta)
    
    # Get baseline sample(s)
    baseline_samples <- rownames(sink_meta[sink_meta$Week == baseline_timepoint, ])
    
    if(length(baseline_samples) > 0) {
      baseline_sample <- baseline_samples[1]  # Take first if multiple
      
      # Calculate distance from baseline to all other timepoints
      for(i in 1:nrow(sink_meta)) {
        current_sample <- rownames(sink_meta)[i]
        
        # Extract distance from distance matrix
        dist_val <- as.matrix(dist_matrix)[baseline_sample, current_sample]
        
        results <- rbind(results, data.frame(
          sink_id = sink_id,
          time_point = sink_meta$Week[i],
          distances = dist_val,
          sample_id = current_sample
        ))
      }
    }
  }
  
  return(results)
}

# Calculate beta diversity from baseline
beta_div_dfs <- beta_div_from_baseline(bray_dist, metadata, baseline_timepoint = 1)

# Save for temporal stability analysis
save(beta_div_dfs, file = "beta_div_graph_dfs.RData")

cat("\n=== Beta Diversity from Baseline Calculated ===\n")
cat("Results saved to: beta_div_graph_dfs.RData\n")

################################################################################
# SAVE RESULTS
################################################################################

# Save distance matrices
saveRDS(aitchison_dist, "aitchison_distance_matrix.rds")
saveRDS(bray_dist, "bray_curtis_distance_matrix.rds")

# Save PERMANOVA results
write.csv(as.data.frame.matrix(permanova_sink), 
          "permanova_overall_aitchison.csv")
write.csv(week_permanova_results, 
          "permanova_by_week_aitchison.csv", 
          row.names = FALSE)

# Save PCoA results
write.csv(pcoa_df, "pcoa_coordinates.csv", row.names = FALSE)

cat("\n=== Beta Diversity Analysis Complete ===\n")
cat("Files saved:\n")
cat("  - aitchison_distance_matrix.rds\n")
cat("  - bray_curtis_distance_matrix.rds\n")
cat("  - permanova_overall_aitchison.csv\n")
cat("  - permanova_by_week_aitchison.csv\n")
cat("  - pcoa_coordinates.csv\n")
cat("  - beta_div_graph_dfs.RData\n")

################################################################################
# NOTES:
# - Aitchison distance: Euclidean distance on CLR-transformed abundances
# - PERMANOVA: Permutational multivariate analysis of variance (999 permutations)
# - CLR transformation handles compositional nature of microbiome data
# - Pseudocount (0.5) added to handle zeros before CLR transformation
# - Bray-Curtis used for temporal stability (distance from baseline timepoint)
################################################################################
