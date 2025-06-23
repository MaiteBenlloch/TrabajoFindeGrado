##########################################################################
# Script Name:    metrics_analysis.R
# Author:         Maite Benlloch Montón
# Date:           [2025-06-23]
# Description:    Computes and compares expression metrics (intensity, variability, 
#                 sparsity) between real and simulated RNA-seq count matrices.
#                 Generates visual comparisons using boxplots and 2D density plots.
#
# Input:          Real and simulated count matrices (genes x samples).
# Output:         Boxplots for expression metrics and 2D density plots of 
#                 intensity vs. sparsity. Global sparsity values are also reported.
#
# Dependencies:   tidyverse, patchwork, scales
#
# Usage Example:  
#   result <- compare_and_plot_metrics(real_counts, simulated_counts, "Dataset Title")
#   print(result$plots)
#   cat("Sparsity - Original:", result$sparsity_global_orig)
#
# Notes:
# - Normalization 
# - Sparsity is defined as the fraction of zero entries.
# - Intensity is the mean normalized count per gene.
# - Variability is computed as the variance across samples.
# - All log transformations use log10(count + 1).
#
##########################################################################


library(patchwork)

# Function to compute metrics from a count matrix
calculate_metrics <- function(counts) {
  # Normalization (Ĉ_ij = C_ij / L_j * median(L_j))
  L_j <- colSums(counts)           # Library sizes
  median_L <- median(L_j)
  counts_norm <- sweep(counts, 2, L_j, "/") * median_L
  
  list(
    intensity = rowMeans(counts_norm),              # Mean expression per gene
    variability = apply(counts_norm, 1, var),       # Variance per gene
    sparsity_global = mean(counts == 0),            # % of zeros in the full matrix
    sparsity_gene = rowMeans(counts == 0),          # % of zeros per gene
    sparsity_sample = colMeans(counts == 0)         # % of zeros per sample
  )
}



spar_glob<- function(original_counts, simulated_counts) {
  stopifnot(all(dim(original_counts) == dim(simulated_counts)))
  
  sparsity_original <- mean(original_counts == 0)
  sparsity_simulated <- mean(simulated_counts == 0)
  
  data.frame(
    Source = c("Original", "Simulated"),
    Sparsity = round(100 * c(sparsity_original, sparsity_simulated), 2)
  )
}

# count_k_ont y rep_k_ont son matrices (o vectores columna) de la muestra "kidney ont"
spar_glob(counts_b_isoseq, rep_b_isoseq)





plot_metrics_from_calculated <- function(counts_original, counts_simulated, sample_name = "Sample") {
  stopifnot(all(dim(counts_original) == dim(counts_simulated)))
  stopifnot(all(colnames(counts_original) == colnames(counts_simulated)))
  
  metrics_orig <- calculate_metrics(counts_original)
  metrics_sim  <- calculate_metrics(counts_simulated)
  
  df_intensity <- data.frame(
    original = metrics_orig$intensity,
    simulated = metrics_sim$intensity
  )
  
  # Para variabilidad, un dataframe largo con etiquetas
  df_variability <- data.frame(
    value = c(metrics_orig$variability, metrics_sim$variability),
    type = rep(c("Original", "Simulated"), each = length(metrics_orig$variability))
  )
  
  # Filtrar outliers extremos para variabilidad (por tipo)
  df_variability_filtered <- do.call(rbind, lapply(split(df_variability, df_variability$type), function(df) {
    Q1 <- quantile(df$value, 0.25)
    Q3 <- quantile(df$value, 0.75)
    IQR <- Q3 - Q1
    df[df$value >= (Q1 - 1.5 * IQR) & df$value <= (Q3 + 1.5 * IQR), ]
  }))
  
  df_sparsity_gene <- data.frame(
    value = c(metrics_orig$sparsity_gene, metrics_sim$sparsity_gene),
    type = rep(c("Original", "Simulated"), each = length(metrics_orig$sparsity_gene))
  )
  
  df_sparsity_sample <- data.frame(
    value = c(metrics_orig$sparsity_sample, metrics_sim$sparsity_sample),
    type = rep(c("Original", "Simulated"), each = length(metrics_orig$sparsity_sample))
  )
  
  p_intensity <- ggplot(df_intensity, aes(x = original, y = simulated)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    coord_cartesian(
      xlim = quantile(df_intensity$original, c(0.01, 0.99)),
      ylim = quantile(df_intensity$simulated, c(0.01, 0.99))
    ) +
    labs(title = "Intensity per Transcript", x = "Original", y = "Simulated") +
    theme_minimal()
  
  p_variability <- ggplot(df_variability_filtered, aes(x = type, y = value, fill = type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    theme_minimal() +
    labs(title = "Variability per Transcript (no extreme outliers)", x = "", y = "Variance") +
    theme(legend.position = "none")
  
  p_sparsity_gene <- ggplot(df_sparsity_gene, aes(x = type, y = value, fill = type)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.size = 0.5) +
    theme_minimal() +
    labs(title = "Sparsity per Transcript", x = "", y = "% Zeros") +
    theme(legend.position = "none")
  
  p_sparsity_sample <- ggplot(df_sparsity_sample, aes(x = type, y = value, fill = type)) +
    geom_boxplot(width = 0.4, outlier.size = 0.5, alpha = 0.7) +
    theme_minimal() +
    labs(title = "Sparsity per Sample", x = "", y = "% Zeros") +
    theme(legend.position = "none")
  
  final_plot <- (p_intensity | p_variability) / (p_sparsity_gene | p_sparsity_sample) +
    plot_annotation(title = paste("Metrics for", sample_name, ": Original vs Simulated"))
  
  return(final_plot)
}



plot_metrics_from_calculated(counts_k_ont, rep_k_ont, sample_name = "Kidney ONT")
plot_metrics_from_calculated(counts_k_isoseq, rep_k_isoseq, sample_name = "Kidney IsoSeq")
plot_metrics_from_calculated(counts_b_ont, rep_b_ont, sample_name = "Brain ONT")
plot_metrics_from_calculated(counts_b_isoseq, rep_b_isoseq, sample_name = "Brain IsoSeq")




#============================================================================
# INTENSITY vs. SPARSITY PLOT
#============================================================================

# Function to calculate intensity and sparsity per gene
intensity_vs_sparsity <- function(counts) {
  L_j <- colSums(counts)
  median_L <- median(L_j)
  counts_norm <- sweep(counts, 2, L_j, "/") * median_L
  
  intensity <- rowMeans(counts_norm)
  sparsity <- rowMeans(counts == 0)
  
  tibble(
    intensity = log10(intensity + 1),  # log10 scale for better visualization
    sparsity = sparsity
  )
}

# Generate data
df_real <- intensity_vs_sparsity(counts_k_isoseq)
df_sim <- intensity_vs_sparsity(rep_k_isoseq)

# Plot: Real Data
p_real <- ggplot(df_real, aes(x = intensity, y = sparsity)) +
  stat_density_2d_filled(contour_var = "ndensity", show.legend = FALSE, geom = "polygon") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(title = "Real Data", x = "Count Intensity (log10)", y = "Sparsity per Transcript [%]") +
  theme_minimal(base_size = 14)

# Plot: Simulated Data
p_sim <- ggplot(df_sim, aes(x = intensity, y = sparsity)) +
  stat_density_2d_filled(contour_var = "ndensity", show.legend = FALSE, geom = "polygon") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(title = "Simulated Data", x = "Count Intensity (log10)", y = "Sparsity per Transcript [%]") +
  theme_minimal(base_size = 14)

# Combine both plots
(p_real | p_sim) +
  plot_annotation(title = "Count Intensity vs. Sparsity per Transcript: Kidney IsoSeq")
#It can be used for the rest of the samples
