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

### Metrics

### Purity testing
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


### Dependency with mean expression testing
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


### Calculation of every clustering metric, but the ROGUE score (ROGUE score is separate)
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


### Calculate ROGUE scores for each cluster for each sample
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


### Correlated gene pairs with Spearman's rho
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