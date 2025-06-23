##########################################################################
# Script Name:    EDA.R
# Author:         Maite Benlloch Montón
# Date:           [2025-06-23]
# Description:    Imports RNA-seq count matrices (ONT and IsoSeq) for brain and 
#                 kidney samples; performs exploratory data analysis including 
#                 zero-count percentages, density plots, and distribution fitting 
#                 to Poisson, Negative Binomial, and Gamma models.
#                 Also analyzes SQANTI structural categories for transcript classification.
#
# Input:          Count matrices in TSV format for ONT and IsoSeq datasets.
# Output:         Plots of zero-count percentages, density distributions, 
#                 distribution fitting results saved as Excel files, and 
#                 barplots of SQANTI structural categories.
#
# Dependencies:   tidyverse, reshape2, ggplot2, writexl, fitdistrplus, dplyr,
#                 patchwork, RColorBrewer, scales
#
# Usage Example:  
#   # Load count matrices and run analysis as shown in the script
#   # Results include saved plots and Excel files with distribution fits
#
# Notes:
# - Count data is rounded and NAs replaced by zero.
# - Distribution fitting is based on AIC to select best model.
# - SQANTI categories are plotted using consistent factor levels.
#
##########################################################################


# IMPORT FILES

#==============================================================================
# KIDNEY - ONT (K31, K32, ... , K35)
#==============================================================================
# Import the TSV file
counts_k_ont <- read.delim("C:/Users/mbenm/Downloads/TFG/ont_data/B0K100/Kconcat/B0K100_counts.tsv", 
                 sep = "\t", 
                 header = TRUE, 
                 stringsAsFactors = FALSE)

# Set the first column as row names
rownames(counts_k_ont) <- counts_k_ont[,1]

counts_k_ont <- counts_k_ont[,-1]
counts_k_ont[is.na(counts_k_ont)] <- 0
counts_k_ont <- round(counts_k_ont)

print("Are columns numeric?")
print(all(sapply(counts_k_ont, is.numeric)))

#==============================================================================
# KIDNEY - ISOSEQ (PACBIO) (K31, K32, ... , K35)
#==============================================================================
counts_k_isoseq <- read.delim("C:/Users/mbenm/Downloads/TFG/isoseq/B0K100/Kconcat/B0K100_counts.tsv", 
                           sep = "\t", 
                           header = TRUE, 
                           stringsAsFactors = FALSE)

# Set the first column as row names
rownames(counts_k_isoseq) <- counts_k_isoseq[,1]
counts_k_isoseq <- counts_k_isoseq[,-1]
counts_k_isoseq[is.na(counts_k_isoseq)] <- 0

counts_k_isoseq <- round(counts_k_isoseq)


print("Are all columns numeric?")
print(all(sapply(counts_k_isoseq, is.numeric)))

#==============================================================================
# BRAIN - ISOSEQ (PACBIO) (B31, B32, ... , B35)
#==============================================================================
counts_b_isoseq <- read.delim("C:/Users/mbenm/Downloads/TFG/isoseq/B100K0/Bconcat/B100K0_counts.tsv", 
                              sep = "\t", 
                              header = TRUE, 
                              stringsAsFactors = FALSE)

# Set the first column as row names
rownames(counts_b_isoseq) <- counts_b_isoseq[,1]
counts_b_isoseq <- counts_b_isoseq[,-1]
counts_b_isoseq[is.na(counts_b_isoseq)] <- 0
counts_b_isoseq <- round(counts_b_isoseq)

print("Are all columns numeric?")
print(all(sapply(counts_b_isoseq, is.numeric)))


#==============================================================================
# BRAIN - ONT (B31, B32, ... , B35)
#==============================================================================
counts_b_ont <- read.delim("C:/Users/mbenm/Downloads/TFG/ont_data/B100K0/Bconcat/B100K0_counts.tsv", 
                              sep = "\t", 
                              header = TRUE, 
                              stringsAsFactors = FALSE)

# Set the first column as row names
rownames(counts_b_ont) <- counts_b_ont[,1]
counts_b_ont <- counts_b_ont[,-1]
counts_b_ont[is.na(counts_b_ont)] <- 0
counts_b_ont <- round(counts_b_ont)

print("Are all columns numeric?")
print(all(sapply(counts_b_ont, is.numeric)))



#==============================================================================
# EXPLORATORY DATA ANALYSIS II
#==============================================================================

# Libraries
library(tidyverse)
library(reshape2)
library(ggplot2)
theme_set(theme_minimal(base_size = 14))  # Base aesthetic

# Function to calculate percentage of zeros per sample
calculate_zero_percentage <- function(count_matrix) {
  zero_percent <- apply(count_matrix, 2, function(col) {
    mean(col == 0) * 100
  })
  data.frame(Sample = names(zero_percent), ZeroPercent = zero_percent)
}

# Function to plot percentage of zeros
plot_zero_percent <- function(zero_df, title) {
  zero_df$Sample <- gsub(".*(K3[1-5]|B3[1-5]).*", "\\1", zero_df$Sample)
  
  ggplot(zero_df, aes(x = Sample, y = ZeroPercent)) +
    geom_bar(stat = "identity", fill = "#4E79A7") +
    theme_gray(base_size = 14) +
    labs(title = title,
         x = "Biological Replicate",
         y = "Percentage of Zero Counts (%)") +
    ylim(0, 100) +
    geom_text(aes(label = sprintf("%.1f%%", ZeroPercent)), vjust = -0.5, size = 8)
}

# Plot density
# For plot_density: transform data before plotting
plot_density <- function(count_matrix, title) {
  df <- as.data.frame(count_matrix)
  df$Transcript <- rownames(count_matrix)
  df_long <- melt(df, id.vars = "Transcript", variable.name = "Sample", value.name = "Count")
  
  df_long$Sample <- gsub(".*(K3[1-5]|B3[1-5]).*", "\\1", df_long$Sample)
  
  df_long$logCount <- log10(df_long$Count + 1)
  
  ggplot(df_long, aes(x = logCount, color = Sample, fill = Sample)) +
    geom_density(alpha = 0.3) +
    theme_gray(base_size = 14) +
    labs(
      title = title,
      x = "log10(Count + 1)",
      y = "Density",
      color = "Biological replicate",
      fill = "Biological replicate"
    )
}


# Set output directory
output_dir <- "C:/Users/mbenm/Downloads/TFG/scripts_def/Visualizations"

# -----------------------------------------------------------------------------
# A. Brain - IsoSeq
zp_b_isoseq <- calculate_zero_percentage(counts_b_isoseq)
p1 <- plot_zero_percent(zp_b_isoseq, "Brain - IsoSeq: Percentage of Zero Counts")
ggsave(filename = file.path(output_dir, "Brain_IsoSeq_ZeroPercent.png"), plot = p1, width = 8, height = 5)

p2 <- plot_density(counts_b_isoseq, "Brain - IsoSeq: Count Distribution (Log10 Scale)")
ggsave(filename = file.path(output_dir, "Brain_IsoSeq_Density.png"), plot = p2, width = 8, height = 5)

# -----------------------------------------------------------------------------
# B. Kidney - IsoSeq
zp_k_isoseq <- calculate_zero_percentage(counts_k_isoseq)
p3 <- plot_zero_percent(zp_k_isoseq, "Kidney - IsoSeq: Percentage of Zero Counts")
ggsave(filename = file.path(output_dir, "Kidney_IsoSeq_ZeroPercent.png"), plot = p3, width = 8, height = 5)

p4 <- plot_density(counts_k_isoseq, "Kidney - IsoSeq: Count Distribution (Log10 Scale)")
ggsave(filename = file.path(output_dir, "Kidney_IsoSeq_Density.png"), plot = p4, width = 8, height = 5)

# -----------------------------------------------------------------------------
# C. Kidney - ONT
zp_k_ont <- calculate_zero_percentage(counts_k_ont)
p5 <- plot_zero_percent(zp_k_ont, "Kidney - ONT: Percentage of Zero Counts")
ggsave(filename = file.path(output_dir, "Kidney_ONT_ZeroPercent.png"), plot = p5, width = 8, height = 5)

p6 <- plot_density(counts_k_ont, "Kidney - ONT: Count Distribution (Log10 Scale)")
ggsave(filename = file.path(output_dir, "Kidney_ONT_Density.png"), plot = p6, width = 8, height = 5)

# -----------------------------------------------------------------------------
# D. Brain - ONT
zp_b_ont <- calculate_zero_percentage(counts_b_ont)
p7 <- plot_zero_percent(zp_b_ont, "Brain - ONT: Percentage of Zero Counts")
ggsave(filename = file.path(output_dir, "Brain_ONT_ZeroPercent.png"), plot = p7, width = 8, height = 5)

p8 <- plot_density(counts_b_ont, "Brain - ONT: Count Distribution (Log10 Scale)")
ggsave(filename = file.path(output_dir, "Brain_ONT_Density.png"), plot = p8, width = 8, height = 5)


#==============================================================================
# EDA III
# Distribution analysis
#==============================================================================


library(writexl)

library(fitdistrplus)

# Function to analyze which distribution fits each sample best
analyze_distributions <- function(count_matrix) {
  results <- data.frame(
    Sample = colnames(count_matrix),
    Best_Distribution = NA,
    AIC_Poisson = NA,
    AIC_NegBin = NA,
    AIC_Gamma = NA
  )
  
  for (i in seq_along(colnames(count_matrix))) {
    sample_name <- colnames(count_matrix)[i]
    data <- count_matrix[, i]
    data <- data[!is.na(data) & data >= 0] # Filter out NA and negative values
    
    # Check if there is enough data to fit distributions
    # At least 2 data points needed to estimate parameters for most distributions.
    # For Gamma distribution, at least 2 data points > 0 are needed for fitdist to work well.
    if (length(data[data > 0]) < 2) { # Only positive values considered for meaningful fit
      results$Best_Distribution[i] <- "Insufficient_Data"
      results$AIC_Poisson[i] <- NA
      results$AIC_NegBin[i] <- NA
      results$AIC_Gamma[i] <- NA
      warning(paste("Sample '", sample_name, "' has insufficient positive data points to fit distributions. Skipping.", sep=""))
      next
    }
    
    # Poisson fit
    # fitdist might throw error if all data are zeros for Poisson.
    # In that case, AIC will be NA and not considered best.
    fit_pois <- tryCatch({
      fitdist(data, "pois")
    }, error = function(e) {
      warning(paste("Error fitting Poisson for sample '", sample_name, "': ", e$message, sep=""))
      NULL
    })
    aic_pois <- if (!is.null(fit_pois)) AIC(fit_pois) else NA
    
    # Negative Binomial fit
    fit_nb <- tryCatch({
      fitdist(data, "nbinom")
    }, error = function(e) {
      warning(paste("Error fitting Negative Binomial for sample '", sample_name, "': ", e$message, sep=""))
      NULL
    })
    aic_nb <- if (!is.null(fit_nb)) AIC(fit_nb) else NA
    
    # Gamma fit (add 1 to avoid zeros, since Gamma distribution requires positive values)
    # fitdist for gamma requires strictly positive values.
    data_gamma <- data + 1
    fit_gamma <- tryCatch({
      fitdist(data_gamma, "gamma")
    }, error = function(e) {
      warning(paste("Error fitting Gamma for sample '", sample_name, "': ", e$message, sep=""))
      NULL
    })
    aic_gamma <- if (!is.null(fit_gamma)) AIC(fit_gamma) else NA
    
    # Select best distribution (lowest AIC)
    aics <- c(Poisson = aic_pois, NegBin = aic_nb, Gamma = aic_gamma)
    # Filter out NA AICs before finding minimum
    valid_aics <- aics[!is.na(aics)]
    
    if (length(valid_aics) > 0) {
      best <- names(which.min(valid_aics))
      results$Best_Distribution[i] <- best
    } else {
      results$Best_Distribution[i] <- "No_Fit_Possible"
    }
    
    results$AIC_Poisson[i] <- aic_pois
    results$AIC_NegBin[i] <- aic_nb
    results$AIC_Gamma[i] <- aic_gamma
  }
  return(results)
}


# Apply the function to each count matrix
results_k_ont     <- analyze_distributions(counts_k_ont)
results_b_ont     <- analyze_distributions(counts_b_ont)
results_k_isoseq  <- analyze_distributions(counts_k_isoseq)
results_b_isoseq  <- analyze_distributions(counts_b_isoseq)

# Print results
print(results_k_ont)
print(results_b_ont)
print(results_k_isoseq)
print(results_b_isoseq)

# Save each table in an Excel file
write_xlsx(results_k_ont,     "C:/Users/mbenm/Downloads/TFG/scripts_def/xlsx/results_kidney_ONT.xlsx")
write_xlsx(results_b_ont,     "C:/Users/mbenm/Downloads/TFG/scripts_def/xlsx/results_brain_ONT.xlsx")
write_xlsx(results_k_isoseq,  "C:/Users/mbenm/Downloads/TFG/scripts_def/xlsx/results_kidney_IsoSeq.xlsx")
write_xlsx(results_b_isoseq,  "C:/Users/mbenm/Downloads/TFG/scripts_def/xlsx/results_brain_IsoSeq.xlsx")



#===============================================================================
# SQANTI Categories
#===============================================================================
library(ggplot2)
library(dplyr)
library(patchwork)
library(RColorBrewer)

files <- c(
  "C:/Users/mbenm/Downloads/TFG/classif_SQANTI/brainisoseq.txt",
  "C:/Users/mbenm/Downloads/TFG/classif_SQANTI/kidneyisoseq.txt",
  "C:/Users/mbenm/Downloads/TFG/classif_SQANTI/brainont.txt",
  "C:/Users/mbenm/Downloads/TFG/classif_SQANTI/kidneyont.txt"
)

# Function to read and count SQANTI categories in each file
read_sqanti_counts <- function(file) {
  data <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)
  table(data$structural_category)
}

# Apply function and combine results
sqanti_counts <- lapply(files, read_sqanti_counts)

# Ensure all tables have same categories (union)
all_categories <- unique(unlist(lapply(sqanti_counts, names)))

sqanti_counts_aligned <- lapply(sqanti_counts, function(tbl) {
  # Add missing categories with zero counts
  missing_cats <- setdiff(all_categories, names(tbl))
  if (length(missing_cats) > 0) {
    tbl <- c(tbl, setNames(rep(0, length(missing_cats)), missing_cats))
  }
  # Order by all_categories
  tbl[all_categories]
})

# Convert to dataframe for plotting
sqanti_df <- do.call(cbind, sqanti_counts_aligned)
colnames(sqanti_df) <- c("Brain IsoSeq", "Kidney IsoSeq", "Brain ONT", "Kidney ONT")
sqanti_df <- as.data.frame(sqanti_df)
sqanti_df$Category <- rownames(sqanti_df)
sqanti_df <- tidyr::pivot_longer(sqanti_df, cols = -Category, names_to = "Sample", values_to = "Count")

# Set factor levels for consistent order
sqanti_df$Category <- factor(sqanti_df$Category, levels = c(
  "full_splice_match",
  "incomplete_splice_match",
  "novel_in_catalog",
  "novel_not_in_catalog",
  "genic",
  "intergenic",
  "fusion",
  "antisense",
  "genic_intron"
))

sqanti_df$Sample <- factor(sqanti_df$Sample, levels = c("Brain IsoSeq", "Kidney IsoSeq", "Brain ONT", "Kidney ONT"))

# Plot SQANTI categories bar plot
ggplot(sqanti_df, aes(x = Sample, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  labs(title = "SQANTI Structural Categories Distribution",
       x = "Sample",
       y = "Number of Transcripts",
       fill = "Category") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set3")
