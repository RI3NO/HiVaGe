# Load required libraries
library(scRNAseq)
library(scater)
library(scDblFinder)
library(scry)
library(SingleR)
library(celldex)

# Load MairPBMCData dataset and store it in 'sce' variable
sce <- MairPBMCData()
sce

# Filtering for empty droplets
# Set random seed for reproducibility
set.seed(42)

# Identify empty droplets using 'emptyDrops' function
e.out <- emptyDrops(counts(sce))

# Summary of empty droplet filtering results
summary(e.out$FDR <= 0.001)

# Create a table showing the counts threshold and FDR filtering results
table(colSums(counts(sce)) > 100, e.out$FDR <= 0.001, useNA = "ifany")

# Apply the filtering based on FDR threshold to 'sce'
sce <- sce[, which(e.out$FDR <= 0.001)]
sce

# Summary of columns and rows with all-zero counts
summary(colMeans(counts(sce) == 0))
summary(rowMeans(counts(sce) == 0))

# Identify cells with all-zero counts for all genes
allzero <- rowMeans(counts(sce) == 0) == 1

# Create a table showing the number of cells with all-zero counts
table(allzero)

# Remove cells with all-zero counts for all genes from 'sce'
sce <- sce[which(!allzero),]

# Log-normalization
# Perform log-normalization on 'sce' using logNormCounts function
sce <- logNormCounts(sce)

# Dimensionality reduction with PCA
# Run PCA on 'sce'
sce <- runPCA(sce)

# Compute doublet density using the reduced dimensions from PCA
dbl.dens <- computeDoubletDensity(sce, d = ncol(reducedDim(sce)))

# Add doublet scores to the 'sce' object
sce$DoubletScore <- dbl.dens

# Summary of doublet scores
summary(dbl.dens)

# Identify doublets using doubletThresholding function
dbl.calls <- doubletThresholding(data.frame(score = dbl.dens),
                                 method = "griffiths", returnType = "call")

# Remove doublets from 'sce' based on the identified doublet calls
sce <- sce[, dbl.calls == "singlet"]

# Perform PCA again
sce <- runPCA(sce)

# GLMPCA clustering
# Perform GLMPCA clustering on 'sce' with L=10 and minibatch="stochastic"
sce <- GLMPCA(sce, L = 10, minibatch = "stochastic")

# tSNE visualization
# Run tSNE on 'sce' using the GLMPCA dimensions
sce <- runTSNE(sce, dimred = "GLMPCA")

# Load reference data for SingleR
ref <- MonacoImmuneData()

# Assign row names to 'sce' using the gene symbols
rownames(sce) <- rowData(sce)$Symbol

# Perform SingleR prediction on 'sce' using the reference data
pred <- SingleR(test = sce, ref = ref, labels = ref$label.fine, assay.type.test = "logcounts")

# Assign the predicted cell types to 'sce$CellTypes'
sce$CellTypes <- factor(pred$pruned.labels)

# Plot tSNE with colors representing cell types
plotTSNE(sce, colour_by = "CellTypes")
