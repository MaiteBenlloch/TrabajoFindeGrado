# ---------------------------------------------------------------
# Script: Simulating Biological Replicates from RNA-seq Count Data
# Description:
# This script prepares RNA-seq count data (counts_k_ont) from a 
# kidney ONT sample, creates a DESeq2 dataset, estimates size 
# factors and dispersions, and simulates biological replicates 
# using Negative Binomial sampling.
# DESeq2 approach
#
# It then compares summary statistics (mean and variance) between 
# the original and simulated counts, visualizing the distributions 
# with overlaid histograms on a log10 scale.
# ---------------------------------------------------------------

#-------------------------------------------------------------
# 1. Load required libraries
#-------------------------------------------------------------
library(DESeq2)

#-------------------------------------------------------------
# 2. Prepare count matrix from counts_k_ont
#-------------------------------------------------------------
counts_k_ont <- as.data.frame(counts_k_ont)  # ensure data.frame

# Set gene names as rownames if needed
if (!is.null(rownames(counts_k_ont))) {
  genes <- rownames(counts_k_ont)
} else {
  genes <- counts_k_ont[[1]]
  counts_k_ont <- counts_k_ont[, -1]
  rownames(counts_k_ont) <- genes
}

# Convert all data to numeric and replace NAs with 0
counts_k_ont[] <- lapply(counts_k_ont, as.numeric)
counts_k_ont[is.na(counts_k_ont)] <- 0  
counts_k_ont <- round(counts_k_ont)

#-------------------------------------------------------------
# 3. Create sample metadata (colData)
#-------------------------------------------------------------
sample_names <- colnames(counts_k_ont)
colData <- data.frame(
  condition = factor(rep("A", length(sample_names))),
  batch     = factor(c("batch1","batch1","batch1","batch2","batch2")),  # adjust if needed
  row.names = sample_names
)

#-------------------------------------------------------------
# 4. Create DESeqDataSet object
#-------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts_k_ont,
                              colData = colData,
                              design  = ~ 1)

#-------------------------------------------------------------
# 5. Estimate size factors and dispersions, then simulate replicates
#-------------------------------------------------------------
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

mu.gene   <- rowMeans(counts(dds, normalized = TRUE))
disp.gene <- dispersions(dds)

valid.genes      <- which(!is.na(mu.gene) & !is.na(disp.gene))
mu.gene.valid    <- mu.gene[valid.genes]
disp.gene.valid  <- disp.gene[valid.genes]
disp.gene.valid[disp.gene.valid == 0] <- 1e-8

set.seed(43)
n.rep <- 5
s.new <- sample(sizeFactors(dds), n.rep, replace = TRUE)

sim.counts <- sapply(s.new, function(sf) {
  rnbinom(n = length(mu.gene.valid),
          mu = mu.gene.valid * sf,
          size = 1 / disp.gene.valid)
})
colnames(sim.counts) <- paste0("sim_rep", seq_len(n.rep))
rownames(sim.counts) <- names(mu.gene.valid)

colData.sim <- data.frame(condition = factor(rep("sim", n.rep)),
                          row.names = colnames(sim.counts))
dds.sim <- DESeqDataSetFromMatrix(countData = sim.counts,
                                  colData = colData.sim,
                                  design = ~ 1)
sizeFactors(dds.sim) <- s.new

#-------------------------------------------------------------
# 6. Prepare real counts for comparison (only simulated genes)
#-------------------------------------------------------------
real_counts <- counts_k_ont[rownames(sim.counts), ]
real_counts <- as.matrix(real_counts)

#-------------------------------------------------------------
# 7. Calculate summary statistics per gene (mean and variance)
#-------------------------------------------------------------
mean_real <- rowMeans(real_counts, na.rm = TRUE)
var_real  <- apply(real_counts, 1, var, na.rm = TRUE)

mean_sim <- rowMeans(sim.counts, na.rm = TRUE)
var_sim  <- apply(sim.counts, 1, var, na.rm = TRUE)

#-------------------------------------------------------------
# 8. Compare distributions with overlaid histograms (log10 scale)
#-------------------------------------------------------------
par(mfrow = c(1, 2))  # split plotting window into 2 panels

# Histogram of means (log10)
hist(log10(mean_real + 1), breaks = 50, col = rgb(0, 0, 1, 0.4),
     main = "Mean distribution (log10)", xlab = "log10(mean + 1)", freq = FALSE)
hist(log10(mean_sim + 1), breaks = 50, col = rgb(1, 0, 0, 0.4),
     add = TRUE, freq = FALSE)
legend("topright", legend = c("Original", "Simulated Biological Replicates"), 
       fill = c(rgb(0,0,1,0.4), rgb(1,0,0,0.4)))

# Histogram of variances (log10)
hist(log10(var_real + 1), breaks = 50, col = rgb(0, 0, 1, 0.4),
     main = "Variance distribution (log10)", xlab = "log10(variance + 1)", freq = FALSE)
hist(log10(var_sim + 1), breaks = 50, col = rgb(1, 0, 0, 0.4),
     add = TRUE, freq = FALSE)
legend("topright", legend = c("Original", "Simulated Biological Replicates"), 
       fill = c(rgb(0,0,1,0.4), rgb(1,0,0,0.4)))

par(mfrow = c(1, 1))  # reset to single panel
