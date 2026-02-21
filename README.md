This repository contains a **comprehensive RNA-seq analysis framework** developed as part of a **Bachelor’s Thesis**. 

The project focuses on:  
- **Exploratory Data Analysis (EDA)** of RNA-seq count matrices 
- **Statistical distribution fitting** (**Poisson**, **Negative Binomial**, **Gamma**) with model selection via **AIC** 
- **Simulation of biological replicates** using multiple modeling strategies (**DESeq2**, **edgeR**, custom mean–variance models) 
- **Sparsity** and **intensity analysis** of gene expression
- **Differential expression simulation** with realistic **biological noise** 
- **Model validation** using **density plots** and **Kolmogorov–Smirnov tests** 

The datasets include **brain** and **kidney** samples sequenced using:  
- **Oxford Nanopore Technologies (ONT)** 
- **PacBio IsoSeq** 

***

## Repository structure

| Script | Description |
|--------|------------|
| **EDA.R** | **Exploratory data analysis** and **distribution fitting** for RNA-seq count matrices. |
| **DESeq2.R** | **Biological replicate simulation** using **DESeq2**’s Negative Binomial model. |
| **MOSim.R** | **Custom simulation** using **mean–variance modeling** (Poisson and Negative Binomial).  |
| **metrics_analysis.R** | Expression metrics comparison (**intensity**, **variability**, **sparsity**). |
| **negativebin_model.R** | Simulation using **edgeR-based Negative Binomial** model with tagwise dispersion. |
| **simulate_DE_noise.R** | **Differential expression simulation** with **biological noise** for benchmarking DE pipelines.  |

***

## Script descriptions

### **EDA.R**

Performs full **exploratory analysis** of RNA-seq count matrices.  

**Features**  
- Imports **ONT** and **IsoSeq** count matrices (brain and kidney).  
- Replaces **NA values** with 0 and rounds counts.  
- Calculates **percentage of zero counts** per replicate (sparsity).  
- Generates **density plots** of log10-transformed counts.  
- Fits **Poisson**, **Negative Binomial**, and **Gamma** distributions to gene counts. 
- Selects the **best-fitting model using AIC**.  
- Exports results as **Excel files**.  
- Visualizes **SQANTI structural transcript categories**.  

**Outputs**  
- Zero-percentage **barplots**.  
- **Density distribution** plots.  
- Excel files with **AIC comparison**.  
- **SQANTI** category stacked barplots.  

***

### **DESeq2.R — DESeq2 approach**

Simulates **biological replicates** using **dispersion** and **size factor** estimates from **DESeq2**’s Negative Binomial model.

**Workflow**  
1. Create a **DESeqDataSet** object.  
2. Estimate **size factors**.  
3. Estimate **dispersion**.  
4. Simulate new replicates via **Negative Binomial sampling**.  

**Comparison**  
- Compares **mean** and **variance** distributions between real and simulated data using **overlaid histograms** (log10 scale).  

**Dependency**  
- **DESeq2 (Bioconductor)**. 

***

### **MOSim.R — Mean–variance modeling approach**

Implements a **custom simulation method** based on a **mean–variance power-law relationship**. 

**Features**  
- Estimates **variance** from mean expression.  
- Uses **Poisson** distribution for low variance and **Negative Binomial** for overdispersed genes. 
- Simulates **biological replicates**.  
- Compares **real vs simulated distributions** via density plots.  

This provides a **simplified parametric alternative** to **DESeq2** and **edgeR** simulation strategies. 

***

### **metrics_analysis.R**

Computes **expression quality metrics** and compares **real vs simulated** datasets.  

**Metrics computed**  
- **Intensity** → mean normalized expression per gene.  
- **Variability** → variance across replicates.  
- **Sparsity (global)** → percentage of zero counts.  
- **Sparsity (per gene)**.  
- **Sparsity (per sample)**.  

**Visualizations**  
- **Scatter plot** (Intensity: original vs simulated).  
- **Boxplots** (variability).  
- **Violin plots** (sparsity per gene).  
- **Boxplots** (sparsity per sample).  
- **2D density plots** (intensity vs sparsity).  

***

### **negativebin_model.R — edgeR approach**

Simulates **biological replicates** using **dispersion estimated from edgeR**.

**Workflow**  
- Estimate **tagwise dispersion** using edgeR’s empirical Bayes method.
- Generate **Negative Binomial replicates**.  
- Adjust **dropout** to match real dataset sparsity.  

**Comparison outputs**  
- **Density plots**.  
- **Per-replicate facet** plots.  
- **Violin plots**.  

**Statistical tests**  
- Performs **Kolmogorov–Smirnov tests**, reporting:  
  - **D statistic**.  
  - **p-value**.  
  - A **similarity conclusion** between real and simulated distributions.  

**Dependency**  
- **edgeR (Bioconductor)**. 

***

### **simulate_DE_noise.R**

Simulates **differential expression** between **control** and **tumor** groups, incorporating **biological noise**.

**Features**  
- Uses **edgeR dispersion estimates**.
- Assigns genes as **upregulated**, **downregulated**, or **non-differential**.  
- Applies predefined **fold changes**.  
- Simulates **biological replicates** under these fold changes.  
- Introduces **Gaussian noise** to model additional biological variability. 

**Visualizations**  
- **Per-replicate density** plots.  
- **Noise vs no-noise** comparisons.  

**Returns**  
- **Simulated counts**.  
- **Noisy counts**.  
- **Fold-change table**.  
- **Density plots**.  

This script enables **realistic benchmarking** of differential expression pipelines under controlled noise conditions. 

***

## Dependencies

**CRAN packages**

```r
install.packages(c(
  "tidyverse", "reshape2", "ggplot2", "dplyr", "tidyr",
  "tibble", "patchwork", "RColorBrewer", "scales",
  "fitdistrplus", "writexl"
))
```

**Bioconductor packages**

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "edgeR"))
```
