# Set-up libraries
library(scRNAseq)
library(DropletUtils)
library(scater)
library(scDblFinder)
library(celldex)
library(SingleR)

library(Rtsne)

library(scry)
library(scran)
library(igraph)

library(clusterSim)

library(fossil)

library(RBGL)

library(tibble)
library(ROGUE)

# Load MairPBMCData dataset and store it in 'sce' variable
sce <- MairPBMCData()
sce

# Create mapping vector of rownames to Symbol
symbol_map = rowData(sce)$Symbol
names(symbol_map) = rownames(sce)

symbol_map

# Fix rownames (duplicates and use of dashes)
rownames(sce) <- rowData(sce)$Symbol
new_rownames = gsub("-", "_", rownames(sce))
rownames(sce) = make.names(new_rownames, unique=TRUE)

standardizeGeneSymbols = function(map, genes) {
  # Voodoo magic that makes all gene names the same across all datasets.. hopefully
  standard_gene_symbols = unname(map[genes])
  standard_gene_symbols[is.na(standard_gene_symbols)] = genes[is.na(standard_gene_symbols)]
  standard_gene_symbols = gsub("-", "_", standard_gene_symbols)  
  return(standard_gene_symbols)
}


### Data preparation (copied from profiling)
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

# Load reference data for SingleR
ref <- MonacoImmuneData()

# Perform SingleR prediction on 'sce' using the reference data
pred <- SingleR(test = sce, ref = ref, labels = ref$label.fine, assay.type.test = "logcounts")

# Assign the predicted cell types to 'sce$CellTypes'
sce$CellTypes <- factor(pred$pruned.labels)



### Metrics
# Retrieve HVGs (M3Drop and seurat)
data <- read.csv("R/flowers_and_clouds/M3Drop_matrix.csv", header = FALSE)
M3Drop_matrix <- as.matrix(data)
M3Drop_matrix <- M3Drop_matrix[2:31, ]
# Standardize names for HVGs
M3Drop_matrix = standardizeGeneSymbols(symbol_map, M3Drop_matrix[,1])

data <- read.csv("R/flowers_and_clouds/seurat_matrix.csv", header = FALSE)
seurat_matrix <- as.matrix(data)
seurat_matrix <- seurat_matrix[2:31, ]
# Standardize names for HVGs
seurat_matrix = standardizeGeneSymbols(symbol_map, seurat_matrix[,1])

### Purity testing
## Define functions
calcPurity = function(clusters, classes) {
  # Just the purity math formula
  sum(apply(table(classes, clusters), 1, max))/length(clusters)
}

getAvrPurity = function(hvgs, sce_object, labels, tSNE_count = 5) {
  print(paste("Chosen t-SNE count:", tSNE_count, sep = ""))
  # Retrieve relevant data
  sce_logcounts = logcounts(sce_object)
  # Vector for storing purity for each clustering
  purities = c(1:tSNE_count)
  
  # Getting expression data for genes of HVGs
  sce_logcounts_genes = as.data.frame(t(sce_logcounts[hvgs,]))
  # Adding cell-type as separate column
  sce_logcounts_genes$CellType = labels
  # Removing duplicates (ignoring the CellType column)
  sce_logcounts_genes = sce_logcounts_genes[!duplicated(subset(sce_logcounts_genes, select=-c(CellType))), ]
  
  for (i in 1:tSNE_count) {
    cat("\r",paste("Running t-SNE #", i, "...", sep = ""))
    # Running t-SNE (with the library from the paper) on expression data
    sce_tSNE = Rtsne(subset(sce_logcounts_genes, select=-c(CellType)), dims = 2) #, perplexity=50, max_iter=2000, early_exag_coeff=12, stop_lying_iter=1000)
    
    # K-means clustering. I chose 28 clusters, since that's the amount of predicted
    # cell types. If that's not how that should be done, pls change
    sce_cluster <- kmeans(sce_tSNE$Y, center=28, nstart=20)
    
    # Calculating purity for (i+1)th clustering
    purities[i] = calcPurity(sce_logcounts_genes$CellType, sce_cluster$cluster)
  }
  
  # Final purity score for this method (seurat) with this database (MairPBMCData)
  purity_score = mean(purities)
  return(list(purities = purities, purity_score = purity_score))
}

## Example use
# Run purity scoring function on HVGs
M3Drop_purity_score = getAvrPurity(M3Drop_matrix, sce, labels = sce$CellTypes)
#seurat_purity_score = getAvrPurity(seurat_matrix, sce, labels = sce$CellTypes)

# First list is 5 purity scores of 5 t-SNE runs, second list is the average value
# and thus it's the final score
M3Drop_purity_score
#seurat_purity_score


### Dependency with mean expression testing
## Define function
getOverlapWithHighLowExpressed <- function(hvgs, sce_object, amount_of_genes_to_check = length(hvgs)) {
  # Get pseudo-bulk data from sce_objectExperiment
  pseudo_bulk = matrix(apply(counts(sce_object), 1, sum))
  rownames(pseudo_bulk) = rownames(counts(sce_object))
  # Log-normalized counts of pseudo-bulk data
  pseudo_bulk = normalizeCounts(pseudo_bulk)
  # Order based on expression
  pseudo_bulk = pseudo_bulk[order(pseudo_bulk, decreasing = TRUE),]
  
  # Get the set amount of highly/lowly expressed genes to test overlap
  hegs = names(pseudo_bulk[1:amount_of_genes_to_check])
  legs = names(pseudo_bulk[(length(pseudo_bulk)-amount_of_genes_to_check):length(pseudo_bulk)])

  # Calculate and return overlap
  overlap_h = table(unique(c(hvgs)) %in% hegs)
  overlap_l = table(unique(c(hvgs)) %in% legs)
  
  expr_var_cor = cor(log(apply(counts(sce_object), 1, mean)), log(apply(counts(sce_object), 1, var)), method = "pearson")
  
  return(list(highly_overlap = unname(overlap_h["TRUE"]/sum(overlap_h)), lowly_overlap = unname(overlap_l["TRUE"]/sum(overlap_l)), correlation_meanVariance = expr_var_cor))
}

#corel = data.frame(var =log(apply(counts(sce), 1, var)), mean = log(apply(counts(sce), 1, mean)))
#ggplot(corel, aes(x = mean, y = var)) + geom_point()+ stat_smooth(method = "lm", col = "red")

## Example use
# Run the dependency test
M3Drop_meanDependence = getOverlapWithHighLowExpressed(M3Drop_matrix, sce)
seurat_meanDependence = getOverlapWithHighLowExpressed(seurat_matrix, sce)

#higly_overlap - overlap with highly expressed genes (positive correlation)
#lowly_overlap - overlap with lowly expressed genes (negative correlation)
M3Drop_meanDependence
seurat_meanDependence


### Calculate CH, DB, AR indexes, Average Silhouette Width and ROGUE score
### based on clustering
## Clustering based on HVGs
set.seed(42)
# Get data from sce_objectExperiment for the HVGs
filtered <- sce[M3Drop_matrix, ]
# Run GLMPCA for dimensionality reduction
filtered <- GLMPCA(filtered, L=10, minibatch="stochastic")
# Build shared nearest-neighbor graph (SNNGraph)
g <- buildSNNGraph(filtered, k=10, use.dimred = 'GLMPCA')
# Perform Louvain clustering on the SNNGraph
clust <- cluster_louvain(g)
# Add cluster data as a factor to "filtered" variable
filtered$Louvain <- factor(membership(clust))
# filtered <- runTSNE(filtered, dimred="GLMPCA")
# plotTSNE(filtered, colour_by="Louvain")

## Define function for calculation of every clustering metric, but the ROGUE score (ROGUE score is separate)
getClusteringMetrics <- function(sce_object, clustering, labels, assay.type = "counts", round = 3){
  # Get relevant data
  if (assay.type == "counts"){
    sce_data = t(counts(sce_object))
  } else if (assay.type == "logcounts"){
    sce_data = t(logcounts(sce_object))
  } else {
    stop("Innappropriate assay type")
  }
  # Make sure clustering and labels are good
  clustering = as.integer(factor(clustering))
  labels = as.integer(factor(labels))
  # Prepare results list
  res = list()
  
  cat("\r","Calculating CH and DB indexes..")
  # Calculate Calinski-Harabasz index
  res$CH_index = round(index.G1(sce_data, clustering), digits=round)
  # Calculate Davies-Bouldin index
  res$DB_index = round(index.DB(sce_data, clustering)$DB, digits=round)
  
  cat("\r","Calculating average silhouette width.. (this one takes awhile)")
  # Calculate average silhouette width
  # First get distance matrix
  sce_dist_matrix = dist(sce_data)
  # Calculate silhouettes
  si_asw = silhouette(clustering, sce_dist_matrix)
  # Get average silhouette width
  res$ASW =  round(summary(si_asw)$avg.width, digits=round)
  
  cat("\r","Calculating Adjusted Rand Index..")
  # Calculate Adjusted Rand Index
  res$ARI_index = round(adj.rand.index(labels, clustering), digits=round)
  
  # Return results
  return(res)
}

## Example use
M3Drop_clustering_metrics = getClusteringMetrics(sce, filtered$Louvain, sce$CellTypes, assay.type = "counts")

M3Drop_clustering_metrics


## Define function to calculate ROGUE scores for each cluster for each sample
getROGUEScore <- function(hvgs, sce_object, clustering, sampling, platform, assay.type = "counts", min.cell.n = 10) {
  # Get all the cluster names
  clusters = unique(clustering)
  # Get all the sample names
  samples = unique(sampling)
  # Prepare results matrix
  res = matrix(nrow = length(clusters), ncol = length(samples))
  rownames(res) = clusters
  colnames(res) = samples
  
  print(paste("Cluster count:", length(clusters), ", sample count:", length(samples), sep = ""))
  
  for (i in 1:length(clusters)) {
    for (j in 1:length(samples)) {
      cat("\r",paste("Working with cluster: ", clusters[i], ", Sample: ", samples[j], sep = ""))
      
      # Get relevant data (cells from cluster i-th and sample j-th)
      if (assay.type == "counts"){
        tmp = counts(sce_object[,clustering == clusters[i] & sampling == samples[j]])
      } else if (assay.type == "logcounts"){
        tmp = logcounts(sce_object[,clustering == clusters[i] & sampling == samples[j]])
      } else {
        stop("Innappropriate assay type")
      }
      
      # Filter out data will too little cells
      if (dim(tmp)[2] >= min.cell.n) {
        # Running the S-E model
        tmp.res <- SE_fun(tmp)
        # Running ROGUE score calculation
        res[i,j] = CalculateRogue(tmp.res, features = hvgs, platform = platform)
      } else {
        # Too little cells!
        res[i,j] = NA
      }
    }
  }
  # Return result matrix
  return(res)
}


## Example use
M3Drop_ROGUE_scores = getROGUEScore(M3Drop_matrix, sce, clustering = filtered$Louvain,
                                    sampling = filtered$Sample_Name, platform = "full-length",
                                    assay.type = "counts")
# Calculating average ROGUE score
M3Drop_avrROGUE = mean(apply(M3Drop_ROGUE_scores, 1, mean, na.rm=TRUE), na.rm=TRUE)
# ROGUE score boxplot
rogue.boxplot(as.data.frame(t(M3Drop_ROGUE_scores)))


### Correlated gene pairs with Spearman's rho
## Define function
correlatedHVGs <- function(hvgs, sce_object) {
  cor_genes <- correlatePairs(sce_object)
  cor_genes <- cor_genes[(cor_genes$gene1 %in% hvgs) & (cor_genes$gene2 %in% hvgs),]
  # Get significant correlated pairs
  sig_cor_genes <- cor_genes$FDR <= 0.05
  #summary(sig_cor_genes)
  # Build undirected graph of correlated genes
  g_cor_genes <- ftM2graphNEL(cbind(cor_genes$gene1, cor_genes$gene2)[sig_cor_genes,],
                              W=NULL, V=NULL, edgemode="undirected")
  # Get clusters of correlated genes
  cl_cor_genes <- highlyConnSG(g_cor_genes)$clusters
  # Order them by size
  cl_cor_genes <- cl_cor_genes[order(lengths(cl_cor_genes), decreasing=TRUE)]
  
  # Prepare result list
  res = list()
  res$gene_correlation_df = cor_genes
  res$correlated_gene_clusters = cl_cor_genes
  # Return results
  return(res)
}

## Example use
M3Drop_correlated_gene_pairs = correlatedHVGs(M3Drop_matrix, sce)

M3Drop_correlated_gene_pairs
