# Temporal Stability and Niche Partitioning of Bacterial Communities in Paired Residential Sink P-traps

This repository contains the code and analytical workflows for the manuscript:

**"Temporal stability and niche partitioning of bacterial communities in paired residential sink P-traps"**

## Overview

This study investigates the microbial communities in paired residential sink P-traps over a 60-day period, examining temporal stability, diversity patterns, and differential abundance of bacterial taxa.

## Repository Structure

```
.
├── 01_preprocessing/          # Data preprocessing and quality control
├── 02_diversity_analysis/     # Alpha and beta diversity analyses
├── 03_differential_abundance/ # ALDEx2 differential abundance analysis
├── 04_temporal_stability/     # Temporal stability and predictive modeling
└── figures/                   # Code for generating publication figures
```

## Data Processing Pipeline

### 1. Preprocessing (Credit: Jorden Rabasco)
- DADA2 amplicon sequence variant (ASV) calling
- Quality filtering and denoising
- Chimera removal
- Taxonomic assignment
- Contamination removal (decontam)
- Rarefaction to 14,700 reads per sample

### 2. Diversity Analysis
- Alpha diversity (Chao1, Shannon indices)
- Beta diversity (Aitchison distance, Bray-Curtis dissimilarity)
- PERMANOVA tests for community composition differences
- Temporal dynamics visualization

### 3. Differential Abundance Analysis
- ALDEx2 analysis with IQLR normalization
- Effect size calculations
- Phenotypic annotation of differentially abundant taxa

### 4. Temporal Stability Analysis
- Coefficient of variation (CV) calculations
- Mean squared successive difference (MSSD)
- Autocorrelation function (ACF) analysis
- Generalized Additive Models (GAMs)
- ARIMA predictive modeling
- Cross-validation for predictive accuracy

## Software Requirements

### R Version
- R 4.3 or higher

### Required R Packages

#### Core packages
```r
dada2 (v1.37.0)
phyloseq (v1.44.0)
vegan (v2.6-4)
ggplot2 (v3.4.4)
dplyr (v1.1.4)
```

#### Diversity and statistics
```r
breakaway
DivNet
mctoolsr
decontam
dunn.test
```

#### Differential abundance
```r
ALDEx2
forcats
```

#### Temporal analysis
```r
mgcv (v1.9)
forecast (v8.21)
```

#### Visualization
```r
ggpubr
tidyr
viridis
cowplot
ggrepel
```

## Installation

Install required packages:

```r
# Install from CRAN
install.packages(c("vegan", "ggplot2", "dplyr", "tidyr", "mgcv", "forecast", 
                   "ggpubr", "viridis", "cowplot", "ggrepel", "forcats", "dunn.test"))

# Install from Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("dada2", "phyloseq", "decontam", "ALDEx2"))

# Install from GitHub
devtools::install_github("adw96/breakaway")
devtools::install_github("adw96/DivNet")
devtools::install_github("leffj/mctoolsr")
```

## Usage

### Running the Complete Analysis

The analyses should be run in sequential order:

```r
# 1. Preprocessing (by Jorden Rabasco)
source("01_preprocessing/dada2_workflow.R")
source("01_preprocessing/contamination_removal.R")

# 2. Diversity analyses
source("02_diversity_analysis/alpha_diversity.R")
source("02_diversity_analysis/beta_diversity.R")

# 3. Differential abundance
source("03_differential_abundance/aldex2_analysis.R")

# 4. Temporal stability
source("04_temporal_stability/stability_metrics.R")
source("04_temporal_stability/predictive_modeling.R")

# 5. Generate figures
source("figures/create_publication_figures.R")
```

## Data Availability

Raw sequencing data will be deposited in the NCBI Sequence Read Archive (SRA) upon publication.

## Citation

If you use this code, please cite:

[Citation information to be added upon publication]

## Contact

For questions about this code repository, please open an issue on GitHub.

## Acknowledgments

- Preprocessing and DADA2 workflow: Jorden Rabasco
- Temporal stability and differential abundance analyses: [Author names]

## License

[License information to be added]
