
# ---------------------------------------------------------------
# Script: Simulating Biological Replicates from RNA-seq Counts
# Description:
# This script simulates biological replicates from RNA-seq count 
# data obtained from a Kidney sample sequenced using Oxford Nanopore 
# Technology (ONT). The simulation is performed using a custom 
# function based on a mean-variance relationship and utilizes 
# Poisson or Negative Binomial distributions as appropriate.
# MOSim approach
#
# The script compares the distribution of the original counts
# with the simulated counts using density plots of log-transformed
# values. This allows visual assessment of how well the simulation
# captures the distributional properties of the real data.
# ---------------------------------------------------------------


# Required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)  # for rownames_to_column

# Define the simulation function
generate_replicates <- function(counts, a = -1, b = 1.2, numReps = 3, minVar = 0.03) {
  counts_clean <- pmax(counts, 0.1, na.rm = TRUE)
  varest <- 10^a * (counts_clean + 1)^b - 1
  varest <- pmax(varest, minVar)
  
  t(apply(cbind(counts_clean, varest), 1, function(x) {
    mean_val <- x[1]
    var_val <- x[2]
    if (mean_val >= var_val) {
      rpois(numReps, lambda = mean_val)
    } else {
      size_param <- mean_val^2 / (var_val - mean_val)
      rnbinom(numReps, size = size_param, mu = mean_val)
    }
  }))
}

# Compute the mean per gene (per row)
gene_means <- rowMeans(counts_k_ont)

# Simulate replicates for each gene with the same number of columns as in the original dataset
set.seed(42)
simulated_matrix <- generate_replicates(gene_means, numReps = ncol(counts_k_ont))
colnames(simulated_matrix) <- paste0("Sim", 1:ncol(counts_k_ont))
rownames(simulated_matrix) <- rownames(counts_k_ont)

# Prepare original data in long format
real_long <- counts_k_ont %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "count") %>%
  mutate(type = "Original")

# Prepare simulated data in long format
sim_long <- as.data.frame(simulated_matrix) %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "count") %>%
  mutate(type = "Simulated biological replicates")

# Combine both datasets
combined <- bind_rows(real_long, sim_long)

# Plot density curves
ggplot(combined, aes(x = log10(count + 1), fill = type, color = type)) +
  geom_density(alpha = 0.4) +
  labs(title = "Kidney ONT: Original vs Simulated (MOSim)",
       x = "log10(Count + 1)",
       y = "Density") +
  theme_minimal()
