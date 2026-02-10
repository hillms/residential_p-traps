################################################################################
# Contamination Removal and Quality Control
# Author: Jorden Rabasco
# Description: Post-DADA2 contamination removal using decontam package and
#              read tracking through the pipeline
################################################################################

# Load required libraries
library(decontam)
library(dplyr)
library(tidyr)
library(ggplot2)

################################################################################
# LOAD DATA
################################################################################

# Load DADA2 objects
pre_dada_stats_list <- readRDS("wills_objs.rds")
dada_results <- readRDS("dada_results.rds")

# Create sequence table
seqtab <- makeSequenceTable(dada_results)

# Remove chimeras
bim2 <- isBimeraDenovo(seqtab, minFoldParentOverAbundance = 4, 
                       multithread = TRUE, verbose = TRUE)

# Report chimera statistics
cat("Percentage of reads identified as chimeric:", 
    sum(seqtab[, bim2])/sum(seqtab) * 100, "%\n")

# Create chimera-free table
temp_bim <- as.data.frame(bim2)
temp_bim <- subset(temp_bim, bim2 != TRUE)
no_chimeras_df <- seqtab[, row.names(temp_bim), drop = FALSE]

################################################################################
# READ METADATA
################################################################################

# Read metadata (AC/BC/FC/PCR indicates control samples)
metadata <- read.csv("m84165_250508_150601_s1.barcodes_summary.csv")

# Identify control samples
metadata$control_or_not <- grepl("AC|BC|PCR|FC", metadata$Sample.Name)

# Identify sink assignment (A or B)
metadata$sink <- ifelse(grepl("_A", metadata$Sample.Name), "A", "B")

# Keep only first 145 samples
metadata <- metadata[1:145, ]
row.names(metadata) <- metadata$Sample.Name

################################################################################
# READ TRACKING THROUGH PIPELINE
################################################################################

# Create read tracking table
read_ret_table <- as.data.frame(
  cbind(
    ccs = pre_dada_stats_list$prim2[, 1],
    primers = pre_dada_stats_list$prim2[, 2], 
    filtered = pre_dada_stats_list$track2[, 2], 
    denoised = sapply(dada_results, function(x) sum(x$denoised))
  )
)

# Add post-chimera counts
read_ret_table$not_chimeras <- rowSums(no_chimeras_df)

# Calculate read retention percentage
read_ret_table$read_ret_per <- read_ret_table$not_chimeras / read_ret_table$ccs

# Add sample names
read_ret_table$Sample_name <- row.names(metadata)

# Remove control samples from tracking
read_ret_table <- read_ret_table[!grepl("C", read_ret_table$Sample_name), ]
row.names(read_ret_table) <- read_ret_table$Sample_name
read_ret_table <- read_ret_table[, 1:5]

################################################################################
# VISUALIZE READ FATE
################################################################################

# Calculate reads removed at each step
graph_df <- read_ret_table
graph_df$rem_primer_removal <- read_ret_table$ccs - read_ret_table$primers
graph_df$rem_filtering <- read_ret_table$primers - read_ret_table$filtered
graph_df$rem_denoised <- read_ret_table$filtered - read_ret_table$denoised
graph_df$rem_chimeras <- read_ret_table$denoised - read_ret_table$not_chimeras
graph_df$reads_retained <- read_ret_table$not_chimeras

# Remove original columns
graph_df <- graph_df[, !colnames(graph_df) %in% 
                       c("ccs", "primers", "filtered", "denoised", 
                         "read_ret_per", "not_chimeras", "Sample_name")]

# Reshape for plotting
sample_names <- rownames(graph_df)
df_long <- graph_df %>%
  pivot_longer(
    cols = colnames(graph_df),
    names_to = "Read_Fate",
    values_to = "Reads"
  )

df_long$sample_names <- unlist(lapply(sample_names, function(x) rep(x, 5)))

# Order categories
ordered_list <- c("rem_primer_removal", "rem_denoised", "rem_filtering",
                  "rem_chimeras", "reads_retained")
df_long <- df_long %>%
  mutate(Read_Fate = factor(Read_Fate, levels = ordered_list))

# Create read fate plot
p_read_fate <- ggplot(df_long, aes(x = sample_names, y = Reads, fill = Read_Fate)) +
  geom_bar(stat = "identity") +
  labs(title = "Read Fate Through Processing Pipeline", 
       x = "Samples", 
       y = "Read Numbers") +
  theme_minimal() +
  theme(axis.text.x = element_blank())

# Save plot
ggsave("read_fate_plot.pdf", p_read_fate, width = 12, height = 6)

# Save read tracking table
write.table(read_ret_table, "read_tracking_table.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE)

################################################################################
# CONTAMINATION REMOVAL WITH DECONTAM
################################################################################

# Run decontam prevalence method
# Threshold = 0.5 for strict contamination removal
decontam_output <- isContaminant(
  no_chimeras_df, 
  neg = metadata$control_or_not, 
  threshold = 0.5, 
  detailed = TRUE, 
  normalize = TRUE, 
  method = 'prevalence'
)

# Plot p-value distribution
pdf("decontam_pvalue_distribution.pdf", width = 8, height = 6)
hist(decontam_output$p, 
     breaks = 20, 
     main = "Distribution of Decontam P-values",
     xlab = "P-value",
     col = "skyblue")
dev.off()

# Report contamination results
n_contaminants <- sum(decontam_output$contaminant)
cat("\nNumber of contaminant ASVs identified:", n_contaminants, "\n")
cat("Percentage of ASVs identified as contaminants:", 
    n_contaminants/nrow(decontam_output) * 100, "%\n")

# Remove contaminant ASVs
true_false_vec <- decontam_output$contaminant != TRUE
no_chimeras_df_clean <- no_chimeras_df[, true_false_vec]

# Report final statistics
cat("\n=== Final ASV Table Statistics ===\n")
cat("Number of samples:", nrow(no_chimeras_df_clean), "\n")
cat("Number of ASVs (after decontam):", ncol(no_chimeras_df_clean), "\n")
cat("Total reads:", sum(no_chimeras_df_clean), "\n")

# Save cleaned ASV table
saveRDS(no_chimeras_df_clean, "asv_table_decontam_clean.rds")
write.table(no_chimeras_df_clean, "asv_table_decontam_clean.txt", 
            sep = "\t", quote = FALSE)

cat("\n=== Contamination Removal Complete ===\n")

################################################################################
# NOTES:
# - Decontam prevalence method compares ASV prevalence in true samples vs controls
# - Threshold = 0.5 is stringent (conservative contamination removal)
# - Control samples identified by patterns: AC, BC, FC, PCR in sample names
# - Read tracking table shows reads retained/removed at each processing step
################################################################################
