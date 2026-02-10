################################################################################
# DADA2 Workflow for Sink Microbiome Study
# Author: Jorden Rabasco
# Description: Amplicon sequence variant (ASV) calling from raw sequencing data
#              using DADA2 pipeline with binned quality error model
################################################################################

# Load required libraries
library(dada2)

# Set working directory (adjust as needed)
# setwd("/path/to/your/data")

################################################################################
# DADA2 PIPELINE
################################################################################

# Step 1: List filtered fastq files
# Files should be primer-trimmed and length-trimmed to 1600bp prior to this step
filts2 <- list.files("/home5/jrabasc/pseudo_pool_dir/data/sink_set/noprimers/filtered", 
                     pattern = "fastq.gz", 
                     full.names = TRUE)

# Step 2: Dereplication
# Combines identical sequences to reduce computational burden
drp2 <- derepFastq(filts2, verbose = TRUE)

# Step 3: Learn Error Rates
# Use binned quality score error model for NovaSeq data
# Binned quality scores at: 3, 10, 17, 22, 27, 35, 40
novaBinnedErrfun <- makeBinnedQualErrfun(c(3, 10, 17, 22, 27, 35, 40))

# Learn error rates from the data
err2 <- learnErrors(drp2, 
                    errorEstimationFunction = novaBinnedErrfun,
                    multithread = TRUE,
                    verbose = TRUE)

# Save error model
saveRDS(err2, file.path("/home5/jrabasc/pseudo_pool_dir/data/sink_set/reg_error_model.rds"))

# Plot error rates for quality control
error_plot <- plotErrors(err2)
saveRDS(error_plot, file.path("/home5/jrabasc/pseudo_pool_dir/data/sink_set/reg_error_model_plot.rds"))

# Step 4: DADA2 Denoising
# Pseudopooling method: processes samples independently but allows rare sequences
# BAND_SIZE=32: controls memory usage and computation speed
dd2 <- dada(drp2, 
            err = err2, 
            pool = "pseudo",
            BAND_SIZE = 32, 
            multithread = TRUE, 
            verbose = TRUE)

# Save denoising results
saveRDS(dd2, file.path("/home5/jrabasc/pseudo_pool_dir/data/sink_set/reg_dada_results.rds"))

# Step 5: Create ASV sequence table
seqtab <- makeSequenceTable(dd2)

# Step 6: Remove chimeras
# minFoldParentOverAbundance=4: stringent chimera detection
bim2 <- isBimeraDenovo(seqtab, 
                       minFoldParentOverAbundance = 4, 
                       multithread = TRUE, 
                       verbose = TRUE)

# Calculate percentage of reads identified as chimeric
chimera_percentage <- sum(seqtab[, bim2]) / sum(seqtab) * 100
cat("Percentage of reads identified as chimeric:", chimera_percentage, "%\n")

# Create data frame of chimera results
temp_bim <- as.data.frame(bim2)
temp_bim <- subset(temp_bim, bim2 != TRUE)

# Create chimera-free sequence table
no_chimeras_df <- seqtab[, row.names(temp_bim), drop = FALSE]

# Save chimera-free ASV table
saveRDS(no_chimeras_df, "asv_table_no_chimeras.rds")

cat("\n=== DADA2 Pipeline Complete ===\n")
cat("Final ASV table dimensions:", dim(no_chimeras_df), "\n")
cat("Number of samples:", nrow(no_chimeras_df), "\n")
cat("Number of ASVs:", ncol(no_chimeras_df), "\n")

################################################################################
# TAXONOMIC ASSIGNMENT (if taxonomy reference database is available)
################################################################################

# Assign taxonomy using SILVA database (version 138.1)
# Download reference from: https://zenodo.org/record/4587955

# taxa <- assignTaxonomy(no_chimeras_df, 
#                        "silva_nr99_v138.1_wSpecies_train_set.fa.gz", 
#                        multithread = TRUE)
# 
# saveRDS(taxa, "taxonomy_assignments.rds")

################################################################################
# NOTES:
# - Raw reads should be primer-trimmed and length-trimmed to 1600bp before this step
# - Binned quality error model is specifically for NovaSeq data
# - Pseudopooling allows detection of rare variants while being computationally efficient
# - minFoldParentOverAbundance=4 provides stringent chimera removal
################################################################################
