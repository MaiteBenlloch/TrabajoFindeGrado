This repository contains a comprehensive RNA-seq analysis framework developed as part of a Bachelor’s Thesis.

The project focuses on:

Exploratory Data Analysis (EDA) of RNA-seq count matrices

Statistical distribution fitting (Poisson, Negative Binomial, Gamma)

Simulation of biological replicates using multiple modeling strategies

Sparsity and intensity analysis

Differential expression simulation with biological noise

Model validation using density plots and Kolmogorov–Smirnov tests

The datasets include brain and kidney samples sequenced using:

Oxford Nanopore Technologies (ONT)

PacBio IsoSeq

Repository Structure
Script	Description
EDA.R	Exploratory data analysis and distribution fitting
DESeq2.R	Biological replicate simulation using DESeq2
MOSim.R	Custom simulation using mean–variance modeling
metrics_analysis.R	Expression metrics comparison (intensity, variability, sparsity)
negativebin_model.R	Simulation using edgeR-based Negative Binomial model
simulate_DE_noise.R	Differential expression simulation with biological noise

Script Descriptions
1.  EDA.R

Performs full exploratory analysis of RNA-seq count matrices.
Features
Imports ONT and IsoSeq count matrices (Brain & Kidney)
Replaces NA values with 0 and rounds counts
Calculates percentage of zero counts per replicate
Generates density plots (log10 scale)
Fits statistical distributions:
Poisson
Negative Binomial
Gamma
Selects best model using AIC
Exports results as Excel files
Visualizes SQANTI structural transcript categories

Outputs:
Zero-percentage barplots
Density distribution plots
Excel files with AIC comparison
SQANTI category stacked barplots

2. DESeq2.R — DESeq2 Approach

Simulates biological replicates using dispersion and size factor estimates from DESeq2.

Workflow

Create DESeqDataSet
Estimate size factors
Estimate dispersion
Simulate new replicates using Negative Binomial sampling

Compare:
Mean distributions
Variance distributions
Visualize results using overlaid histograms (log10 scale)

Dependency
DESeq2 (Bioconductor)

3. MOSim.R — Mean–Variance Modeling Approach

Custom simulation method based on a mean–variance power-law relationship.

Features
Estimates variance from mean expression
Uses:
Poisson distribution (low variance)
Negative Binomial distribution (overdispersion)
Simulates biological replicates
Compares real vs simulated distributions via density plots
This approach provides a simplified parametric alternative to DESeq2 and edgeR.

4. metrics_analysis.R

Computes expression quality metrics and compares real vs simulated datasets.

Metrics Computed
Intensity → Mean normalized expression per gene
Variability → Variance across replicates
Sparsity (global) → % of zero counts
Sparsity (per gene)
Sparsity (per sample)

Visualizations
Scatter plot (Intensity Original vs Simulated)
Boxplots (Variability)
Violin plots (Sparsity per gene)
Boxplots (Sparsity per sample)
2D density plots (Intensity vs Sparsity)

5. negativebin_model.R — edgeR Approach

Simulates biological replicates using dispersion estimated from edgeR.

Workflow
Estimate tagwise dispersion (edgeR)
Generate Negative Binomial replicates
Adjust dropout to match real sparsity

Compare distributions:
Density plots
Per-replicate facet plots
Violin plots
Perform Kolmogorov–Smirnov test:
Reports D statistic
Reports p-value
Provides similarity conclusion

Dependency: edgeR

6. simulate_DE_noise.R

Simulates differential expression between control and tumor groups, incorporating biological noise.

Features
Uses edgeR dispersion estimates
Assigns genes as:
Upregulated
Downregulated
Non-differential
Applies fold changes
Simulates biological replicates
Introduces Gaussian noise

Visualizes:
Per-replicate density plots
Noise vs no-noise comparison

Returns:
Simulated counts
Noisy counts
Fold change table
Density plots
This script allows realistic benchmarking of DE pipelines.

Dependencies

Install required CRAN packages:
install.packages(c(
  "tidyverse",
  "reshape2",
  "ggplot2",
  "dplyr",
  "tidyr",
  "tibble",
  "patchwork",
  "RColorBrewer",
  "scales",
  "fitdistrplus",
  "writexl"
))

Install Bioconductor packages:

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "edgeR"))
