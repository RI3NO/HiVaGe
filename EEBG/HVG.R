#Lib loading
# Load necessary libraries
library(scRNAseq)
library(scater)
library(scDblFinder)
library(scry)
library(BASiCS)
library(igraph)
library(GGally)
library(scran)
library(remotes)
library(Seurat)
library(tibble)
library(ROGUE)
library(DropletUtils)
library(M3Drop)
library(Seurat)
library(SingleR)
library(celldex)
install_github("hillas/scVEGs")

# Load MairPBMCData dataset and store it in 'sce' variable
sce <- MairPBMCData()
sce

# Filtering
# Remove empty droplets from the dataset
set.seed(42)
e.out <- emptyDrops(counts(sce))
summary(e.out$FDR <= 0.001)

# Table summarizing the filtering results
table(colSums(counts(sce)) > 100, e.out$FDR <= 0.001, useNA = "ifany")

# Apply the filtering to 'sce'
sce <- sce[, which(e.out$FDR <= 0.001)]
sce

# Check for cells with all-zero counts and remove them
summary(colMeans(counts(sce) == 0))
summary(rowMeans(counts(sce) == 0))
allzero <- rowMeans(counts(sce) == 0) == 1
table(allzero)
sce <- sce[which(!allzero),]

# Log normalization
# Perform log-normalization on 'sce'
sce <- logNormCounts(sce)

# Run PCA on 'sce'
sce <- runPCA(sce)

# Compute doublet density
dbl.dens <- computeDoubletDensity(sce, d = ncol(reducedDim(sce)))
sce$DoubletScore <- dbl.dens
summary(dbl.dens)

# Identify and remove doublets based on thresholding
dbl.calls <- doubletThresholding(data.frame(score = dbl.dens),
                                 method = "griffiths", returnType = "call")
sce <- sce[, dbl.calls == "singlet"]

# Clustering
# Run GLMPCA for dimensionality reduction
set.seed(42)
sce <- GLMPCA(sce, L = 10, minibatch = "stochastic")

# Build shared nearest-neighbor graph (SNNGraph)
g <- buildSNNGraph(sce, k = 10, use.dimred = 'GLMPCA')
g

# Perform Louvain clustering on the SNNGraph
clust <- igraph::cluster_louvain(g)
sce$Louvain <- factor(membership(clust))

# BASiCS
# Create BASiCS data object 'DataNoSpikes'
DataNoSpikes <- newBASiCS_Data(
  counts(sce), Tech = FALSE, SpikeInfo = NULL,
  BatchInfo = c14$Sample_Name
)
DataNoSpikes

# Run BASiCS MCMC
ChainNoSpikes <- BASiCS_MCMC(
  Data = DataNoSpikes, N = 1000,
  Thin = 10, Burn = 500,
  WithSpikes = FALSE, Regression = TRUE,
  PrintProgress = FALSE
)

# Detect highly variable genes (HVGs) using BASiCS
HVG <- BASiCS_DetectVG(
  ChainNoSpikes,
  Task = c("HVG"),
  PercentileThreshold = 0.2,
  VarThreshold = NULL,
  ProbThreshold = 0.5,
  EpsilonThreshold = NULL,
  EFDR = 0.1,
  Plot = FALSE,
  MinESS = 100,
)

# Store the HVGs in 'BASiCS_matrix'
BASiCS_matrix <- matrix(nrow = 30, ncol = 13)
BASiCS_genes <- unlist(as.data.frame(HVG)$GeneName)
BASiCS_matrix[, 1] <- BASiCS_genes
BASiCS_matrix

# Loop through clusters and run BASiCS for each cluster
cluster_list <- list("c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13")
BASiCS_matrix <- matrix(nrow = 30, ncol = length(cluster_list))

for (i in 1:length(cluster_list)) {
  # Get the dataset corresponding to the cluster
  dataset_name <- cluster_list[[i]]
  dataset <- get(dataset_name)
  # Run BASiCS on the dataset
  # ... (omitted for brevity)
  # Store the result in the BASiCS_matrix
  # ... (omitted for brevity)
}

# Export the BASiCS_matrix to CSV
write.csv(BASiCS_matrix, file = 'BASiCS_matrix.csv', row.names = FALSE)

# Scran export
# Generate a new Seurat object 'dec.pbmc' for modeling gene variability
set.seed(42)
dec.pbmc <- modelGeneVar(sce)

# Select the top highly variable genes (HVGs) using scran
scran_genes <- getTopHVGs(sce, n = 200)

# Create a matrix 'scran_matrix' to store the top HVGs
scran_matrix <- matrix(nrow = 30, ncol = 13)

# Store the scran genes in the first column of the scran_matrix
scran_genes <- unlist(top.pbmc)[1:30]
scran_matrix[, 1] <- scran_genes
scran_matrix

# Loop scran
# Iterate through the clusters and perform scran on each cluster
for (i in 1:14) {
  # Generate the variable name based on the iteration
  var_name <- paste0("c", i)
  
  # Execute the code snippet with the specific variable name
  eval(parse(text = paste0("
    dec.pbmc <- modelGeneVar(", var_name, ")
    top.pbmc <- getTopHVGs(", var_name, ", n=30)
    scran_genes <- unlist(top.pbmc)[1:30]
    scran_matrix[, ", i, "] <- scran_genes
  ")))
}

# Print the resulting matrix
scran_matrix

# Scran export
# Export the scran_matrix to a CSV file
write.csv(scran_matrix, file = 'scran_matrix.csv', row.names = FALSE)


# M3Drop
# Identify highly variable genes (HVGs) using M3Drop
M3Drop_genes <- counts(sce) %>%
  NBumiConvertData() %>%
  NBumiFitModel() %>%
  NBumiFeatureSelectionCombinedDrop(ntop = 200) %>%
  rownames()
M3Drop_genes

# Create a matrix 'M3Drop_matrix' to store the M3Drop HVGs
M3Drop_matrix <- matrix(nrow = 30, ncol = 13)

# Store the M3Drop genes in the first column of the M3Drop_matrix
M3Drop_genes <- unlist(variable_genes_M3Drop)[1:30]
M3Drop_matrix[, 1] <- M3Drop_genes
M3Drop_matrix


# ROGUE
# Run ROGUE algorithm to detect differentially expressed genes
ROGUE <- SE_fun(count_matrix, span = 0.5, r = 1, mt.method = "fdr", if.adj = T)
ROGUE$Gene

# Create a matrix 'ROGUE_matrix' to store the ROGUE DE genes
ROGUE_matrix <- matrix(nrow = 30, ncol = 13)

# Store the ROGUE genes in the first column of the ROGUE_matrix
ROGUE_genes <- unlist(ROGUE$Gene)[1:200]
ROGUE_matrix[, 1] <- ROGUE_genes
ROGUE_matrix

# ROGUE loop
# Iterate through the clusters and run ROGUE on each cluster
ROGUE_matrix <- matrix(nrow = 30, ncol = length(cluster_list))

for (i in 1:length(cluster_list)){
  dataset_name <- cluster_list[[i]]
  dataset <- get(dataset_name)
  count_matrix <- assay(c1, "counts")
  ROGUE <- SE_fun(count_matrix, span = 0.5, r = 1, mt.method = "fdr", if.adj = T)
  ROGUE$Gene
  # Store the result in the ROGUE_matrix
  ROGUE_genes <- unlist(ROGUE$Gene)[1:30]
  ROGUE_matrix[, i] <- ROGUE_genes
}

# Export the ROGUE_matrix to a CSV file
write.csv(ROGUE_matrix, file = 'ROGUE_matrix.csv', row.names = FALSE)


# Extract the count matrix from the SCE object
count_matrix <- assay(sce, "counts")

# Convert the count matrix to a regular matrix
count_matrix <- as.matrix(count_matrix)
rownames(count_matrix) <- rownames(sce)
colnames(count_matrix) <- colnames(sce)

# loop M3Drop
# Iterate through the clusters and perform M3Drop on each cluster
cluster_list <- list("c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12", "c13")
M3Drop_matrix <- matrix(nrow = 30, ncol = length(cluster_list))

for (i in 1:length(cluster_list)) {
  dataset_name <- cluster_list[[i]]
  dataset <- get(dataset_name)
  variable_genes_M3Drop <- counts(sce) %>%
    NBumiConvertData() %>%
    NBumiFitModel() %>%
    NBumiFeatureSelectionCombinedDrop(ntop = 30) %>%
    rownames()
  
  variable_genes_M3Drop
  # Store the result in the M3Drop_matrix
  M3Drop_genes <- variable_genes_M3Drop
  M3Drop_matrix[, i] <- M3Drop_genes
}

# Print the resulting matrix
M3Drop_matrix

# M3Drop export
# Export the M3Drop_matrix to a CSV file
write.csv(M3Drop_matrix, file = 'M3Drop_matrix.csv', row.names = FALSE)


# Seurat
# Create a data frame 'data' from the log-normalized counts of 'sce'
data <- as.data.frame(logcounts(sce))

# Create a Seurat object 'seurat_obj' from the data frame
seurat_obj <- CreateSeuratObject(data)

# Scale the Seurat object
seurat_obj <- ScaleData(seurat_obj)

# Find variable features using Seurat's VST method
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 200, verbose = FALSE)

# Get the variable genes identified by Seurat
seurat_genes <- seurat_obj@assays$RNA@var.features
seurat_genes

# Viz tSNE clusters
# Assign the Louvain clustering results to 'sce$Louvain'
sce$Louvain <- colData$Louvain

# Perform tSNE analysis on 'sce' using the GLMPCA dimensions
sce <- runTSNE(sce, dimred = "GLMPCA")

# Load reference data for SingleR
ref <- MonacoImmuneData()

# Perform SingleR prediction on 'sce' using the reference data
pred <- SingleR(test = sce, ref = ref, labels = ref$label.fine, assay.type.test = "logcounts")

# Assign the predicted cell types to 'sce$CellTypes'
sce$CellTypes <- factor(pred$pruned.labels)

# Plot tSNE with colors representing cell types
plotTSNE(sce, colour_by = "CellTypes")

# ROGUE_genes, BASiCS_genes, seurat_genes, scran_genes, M3Drop_genes

# Create a Venn diagram to visualize the overlap of whole gene sets obtained by different methods, not taking into account clusters
venn.diagram(
  x = list(BASiCS_genes, M3Drop_genes, ROGUE_genes, scran_genes, seurat_genes),
  category.names = c("BASiCS", "M3Drop", "ROGUE", "scran", "seurat"),
  filename = "Whole_Venn_diagram.png",
  output = TRUE,
  col = color_palette,   # Apply custom color palette
  fill = color_palette,  # Fill circles with custom colors
  alpha = 0.5,           # Adjust transparency of circles
  cex = 1.5              # Increase font size of labels
)

# Clustering based on HVGs

# Use runPCA to perform PCA on 'sce' again
sce <- runPCA(sce)

# Plot PCA to visualize clustering patterns
plotPCA(sce)

# Filter 'sce' using M3Drop HVGs
filtered <- sce[M3Drop_genes, ]

# Perform GLMPCA on the filtered dataset
filtered <- GLMPCA(filtered, L = 10, minibatch = "stochastic")

# Run tSNE on the filtered dataset using the GLMPCA dimensions
sce <- runTSNE(sce, dimred = "GLMPCA")

# Load reference data for SingleR
ref <- MonacoImmuneData()

# Assign row names to 'sce'
rownames(sce) <- rowData(sce)$Symbol

# Perform SingleR prediction on the filtered dataset using the reference data
pred_filtered <- SingleR(test = filtered, ref = ref, labels = ref$label.fine, assay.type.test = "logcounts")

# Assign the predicted cell types to 'filtered$CellTypes'
filtered$CellTypes <- factor(pred_filtered$pruned.labels)

# Plot tSNE with colors representing cell types for the filtered dataset
plotTSNE(filtered, colour_by = "CellTypes")

# Calculate the percentage of cells for each cell type in the whole dataset and filtered dataset
levels <- unique(sce$CellTypes)
cell_percentages <- matrix(1:2 * length(levels), nrow = 2, ncol = length(levels), byrow = TRUE)
levels <- as.vector(levels)
colnames(cell_percentages) <- levels
rownames(cell_percentages) <- c("whole", "HVGs") 
filtered$CellTypes

# Calculate cell percentages for the whole dataset
for (i in 1:length(levels)) {
  count <- table(sce$CellTypes == levels[i])["TRUE"]
  cell_percentages[1, i] <- count / length(sce$CellTypes) * 100
}

# Calculate cell percentages for the filtered dataset
for (i in 1:length(levels)) {
  count <- table(filtered$CellTypes == levels[i])["TRUE"]
  cell_percentages[2, i] <- count / length(filtered$CellTypes) * 100
}

# Export cell_percentages to a CSV file
write.csv(cell_percentages, file = "cell_percentages.csv", row.names = TRUE)
