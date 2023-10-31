calculate_RDR <- function(genes_list1, genes_list2) {
  intersect_genes <- intersect(genes_list1, genes_list2)
  RDR <- length(intersect_genes) / length(genes_list1)
  return(RDR)
}

setwd("C:/Users/redegator/PycharmProjects/HiVaGe")
#Lib loading
# Load necessary libraries
library(scRNAseq)
library(scater)
library(scDblFinder)
library(scran)
library(scry)
library(BASiCS)
library(DropletUtils)
library(Seurat)

# Load MairPBMCData dataset and store it in 'sce' variable
sce <- MairPBMCData()
sce

# Set the random seed for reproducibility
set.seed(42)

# Number of cells in the dataset
n_cells <- nrow(sce)

# Randomly shuffle the row indices to create a random order
random_order <- sample(n_cells)

# Calculate the midpoint to split the dataset into two equal groups
midpoint <- n_cells / 2

# Divide the shuffled indices into two groups (training and validation)
training_indices <- random_order[1:midpoint]
validation_indices <- random_order[(midpoint + 1):n_cells]

# Create the training and validation datasets
sce_training <- sce[training_indices, ]
sce_validation <- sce[validation_indices, ]

# Filtering
# Remove empty droplets from the dataset
set.seed(43)
e.out_training <- emptyDrops(counts(sce_training))
e.out_validation <- emptyDrops(counts(sce_validation))
summary(e.out_training$FDR <= 0.001)
summary(e.out_validation$FDR <= 0.001)

# Table summarizing the filtering results
table(colSums(counts(sce_training)) > 100, e.out_training$FDR <= 0.001, useNA = "ifany")
table(colSums(counts(sce_validation)) > 100, e.out_validation$FDR <= 0.001, useNA = "ifany")

# Apply the filtering to 'sce'
sce_training <- sce[, which(e.out_training$FDR <= 0.001)]
sce_validation <- sce[, which(e.out_validation$FDR <= 0.001)]
sce_training
sce_validation

# Check for cells with all-zero counts and remove them
summary(colMeans(counts(sce_training) == 0))
summary(rowMeans(counts(sce_training) == 0))
allzero_training <- rowMeans(counts(sce_training) == 0) == 1
table(allzero_training)
sce_training <- sce_training[which(!allzero_training),]

summary(colMeans(counts(sce_validation) == 0))
summary(rowMeans(counts(sce_validation) == 0))
allzero_validation <- rowMeans(counts(sce_validation) == 0) == 1
table(allzero_validation)
sce_validation <- sce_validation[which(!allzero_validation),]

# Log normalization
# Perform log-normalization on 'sce'
sce_training <- logNormCounts(sce_training)
sce_validation <- logNormCounts(sce_validation)

# Run PCA on 'sce'
sce_training <- runPCA(sce_training)
sce_validation <- runPCA(sce_validation)

# Compute doublet density
dbl.dens_training <- computeDoubletDensity(sce_training, d = ncol(reducedDim(sce_training)))
sce_training$DoubletScore <- dbl.dens_training
summary(dbl.dens_training)

dbl.dens_validation <- computeDoubletDensity(sce_validation, d = ncol(reducedDim(sce_validation)))
sce_validation$DoubletScore <- dbl.dens_validation
summary(dbl.dens_validation)

# Identify and remove doublets based on thresholding
dbl.calls_training <- doubletThresholding(data.frame(score = dbl.dens_training),
                                 method = "griffiths", returnType = "call")
sce_training <- sce_training[, dbl.calls_training == "singlet"]

dbl.calls_validation <- doubletThresholding(data.frame(score = dbl.dens_validation),
                                          method = "griffiths", returnType = "call")
sce_validation <- sce_validation[, dbl.calls_validation == "singlet"]
sce_training
sce_validation


# Scran 
# Generate a new Seurat object 'dec.pbmc' for modeling gene variability
set.seed(43)
dec.pbmc_training <- modelGeneVar(sce_training)
dec.pbmc_validation <- modelGeneVar(sce_validation)

# Select the top highly variable genes (HVGs) using scran
scran_genes_training <- getTopHVGs(dec.pbmc_training, n = 200)
scran_genes_validation <- getTopHVGs(dec.pbmc_validation, n = 200)
# scran_genes_training
# scran_genes_validation

# Calculate RDR
calculate_RDR(scran_genes_training, scran_genes_validation)


# Seurat
# Create a data frame 'data' from the log-normalized counts of 'sce'
data_training <- as.data.frame(logcounts(sce_training))
data_validation <- as.data.frame(logcounts(sce_validation))

# Create a Seurat object 'seurat_obj' from the data frame
seurat_obj_training <- CreateSeuratObject(data_training)
seurat_obj_validation <- CreateSeuratObject(data_validation)

# Scale the Seurat object
seurat_obj_training <- ScaleData(seurat_obj_training)
seurat_obj_validation <- ScaleData(seurat_obj_validation)

# Find variable features using Seurat's VST method
seurat_obj_training <- FindVariableFeatures(seurat_obj_training, selection.method = "vst", nfeatures = 200, verbose = FALSE)
seurat_obj_validation <- FindVariableFeatures(seurat_obj_validation, selection.method = "vst", nfeatures = 200, verbose = FALSE)

# Get the variable genes identified by Seurat
seurat_genes_training <- seurat_obj_training@assays$RNA@var.features
seurat_genes_validation <- seurat_obj_validation@assays$RNA@var.features
# seurat_genes_training
# seurat_genes_validation

# Calculate RDR
calculate_RDR(seurat_genes_training, seurat_genes_validation)

