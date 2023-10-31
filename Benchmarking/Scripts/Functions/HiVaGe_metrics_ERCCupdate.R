### Set-up libraries for metrics
library(dplyr)
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

### Set-up libraries for HVGs
library(reticulate)
source_python("HiVaGePY.py")

library(magrittr)
library(Seurat)
source("scVEGs.R")
library(singleCellHaystack)
library(scmap)
library(SIEVE)
library(scLVM)
library(tictoc)

## Message printing function for outputting progress
replaceCat <- function(x, width = 50)
{
  cat("\r",rep(" ", times = width - length(x)), "\r")
  cat(x)
}


### Testing the dependency with the mean expression
# hvgs - vector containing HVGs
# sce_object - SigleCellExperiment object
# batching=1 - vector containing batches' number for each cell. doesn't impact the result if not set
# amount_of_genes_to_check=length(hvgs) - the amount of highly expressed genes to check the overlap of HVGs with. if not set, equals the length of hvgs vector
getOverlapWithHighLowExpressed <- function(hvgs, sce_object, batching = 1, amount_of_genes_to_check = length(hvgs)) {
  # Processing batching input
  batches = levels(factor(batching))
  # Calculation of pseudobulk expression data
  pseudo_bulk = matrix(nrow = length(rownames(sce_object)), ncol = length(batches))
  colnames(pseudo_bulk) = batches
  for (i in batches) {
    pseudo_bulk[,i] = apply(counts(sce_object[, batching == i]), 1, sum)
  }
  # Normalization between different batches
  pseudo_bulk = normalizeCounts(pseudo_bulk)
  # Vector of average explassion of the genes
  avr_expr = apply(pseudo_bulk, 1, mean)
  names(avr_expr) = rownames(counts(sce_object))
  # Sorting vector's values from highest expression to lowest
  avr_expr = sort(avr_expr, decreasing = TRUE)
  
  # Get the set amount of highly/lowly expressed genes to test overlap
  hegs = names(avr_expr[1:amount_of_genes_to_check])
  legs = names(avr_expr[(length(avr_expr)-amount_of_genes_to_check):length(avr_expr)])
  
  # Calculate and return overlap
  overlap_h = table(unique(c(hvgs)) %in% hegs)
  overlap_l = table(unique(c(hvgs)) %in% legs)

  # Calculate general correlation between gene expression varience and the mean
  expr_var_cor = cor(log(apply(counts(sce_object), 1, mean)), log(apply(counts(sce_object), 1, var)), method = "pearson")

  # Return list containing overlap with highly expressed genes, lowly expressed genes and correaltion mean~variance for the dataset.
  return(list(highly_overlap = unname(overlap_h["TRUE"]/sum(overlap_h)), lowly_overlap = unname(overlap_l["TRUE"]/sum(overlap_l)), correlation_meanVariance = expr_var_cor))
}


### Purity testing
## Function for calculatioon of purity in accordance with its math formula 
# clusters - vector containing cluster's number for each cell
# classes - vector containing cell type for each cell
calcPurity = function(clusters, classes) {
  sum(apply(table(classes, clusters), 1, max))/length(clusters)
}

## Function for calculating average purity
getAvrPurity = function(hvgs, sce_object, labels, tSNE_count = 5) {
  # Setting seed for consistent clustering
  set.seed(42)
  # Outputting information about the running function
  replaceCat(paste("Chosen t-SNE count:", tSNE_count, ".\n", sep = ""))
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

  # Disallowing for division by 0 further down the line
  if (nrow(sce_logcounts_genes) > 3) {
    
    # Fix the perplexity with accordance to amount of genes
    moving_perplexity = floor((nrow(sce_logcounts_genes) - 2) / 3)
    if (moving_perplexity > 50) {
      moving_perplexity = 50
    }

    # Repeat tSNE the required set of times
    for (i in 1:tSNE_count) {
      replaceCat(paste("Running t-SNE #", i, "...", sep = ""))
      # Running t-SNE (with the library from the paper) on expression data
      sce_tSNE = Rtsne(subset(sce_logcounts_genes, select=-c(CellType)), dims = 2, check_duplicates = FALSE, perplexity=moving_perplexity) # max_iter=2000, early_exag_coeff=12, stop_lying_iter=1000)
      
      cluster_num = length(levels(factor(labels)))
      
      # K-means clustering. I chose 28 clusters, since that's the amount of predicted
      # cell types. If that's not how that should be done, pls change
      sce_cluster <- kmeans(sce_tSNE$Y, center=cluster_num, nstart=20)
      
      # Calculating purity for (i+1)th clustering
      purities[i] = calcPurity(sce_logcounts_genes$CellType, sce_cluster$cluster)
    }
    
    # Final purity score for this method (seurat) with this database (MairPBMCData)
    purity_score = mean(purities)
    return(purity_score)
  } else {
    replaceCat(paste("Not enough HVGs, skipping purity score.\n"))
  }
}


### Calculate ROGUE scores for each cluster for each sample
# hvgs - vector containing HVGs
# sce_object - SigleCellExperiment object
# clustering - vector containing cluster's number for each cell
# sampling - vector containing sample's number for each cell. doesn't impact the result if not set
# platform - the platform ("UMI" or "full-length") used for generating the tested dataset
# assay.type - either "counts" or "logcounts"
# min.cell.n = 10 - minimum cell count allowed for ROGUE score calculation
getROGUEScore <- function(hvgs, sce_object, clustering, sampling, platform, assay.type = "counts", min.cell.n = 10) {
  # Get all the cluster names
  clusters = unique(clustering)
  # Get all the sample names
  samples = unique(sampling)
  # Prepare results matrix
  res = matrix(nrow = length(clusters), ncol = length(samples))
  rownames(res) = clusters
  colnames(res) = samples

  # Print info about the running function
  replaceCat(paste("Cluster count:", length(clusters), ", sample count:", length(samples), sep = ""))
  
  for (i in 1:length(clusters)) {
    for (j in 1:length(samples)) {
      # Print progress
      replaceCat(paste("Working with cluster: ", clusters[i], ", Sample: ", samples[j], sep = ""))
      
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


### Calculate correlation between HVGs
# hvgs - vector containing HVGs
# filtered_sce_object - SingleCellExperiment, that only has HVGs in it
correlatedHVGs <- function(hvgs, filtered_sce_object) {
  # Calculate correlation between paris of HVGs
  cor_genes <- correlatePairs(filtered_sce_object)
  cor_genes <- cor_genes[!is.na(cor_genes$FDR),]
  # Prepare result list
  res = list()
  res$gene_correlation_df = cor_genes
  res$avr_abs_rho = mean(abs(cor_genes$rho))

  #cor_genes <- cor_genes[(cor_genes$gene1 %in% hvgs) & (cor_genes$gene2 %in% hvgs),]
  # Get significant correlated pairs
  sig_cor_genes <- cor_genes$FDR <= 0.05
  #summary(sig_cor_genes)
  # Build undirected graph of correlated genes
  if (!is.null(dim(cbind(cor_genes$gene1, cor_genes$gene2)[sig_cor_genes,]))) {
    g_cor_genes <- ftM2graphNEL(cbind(cor_genes$gene1, cor_genes$gene2)[sig_cor_genes,],
                                W=NULL, V=NULL, edgemode="undirected")
    # Get clusters of correlated genes
    cl_cor_genes <- highlyConnSG(g_cor_genes)$clusters
    # Order them by size
    cl_cor_genes <- cl_cor_genes[order(lengths(cl_cor_genes), decreasing=TRUE)]
    
    res$correlated_gene_clusters = cl_cor_genes
  }
  
  # Return results
  return(res)
}


### Get HVGs based on method
HiVaGe <- function(sce_object, flavour, batching = 1, num_HVGs = 400, ERCCs = NULL) {
  # Should already be preinstalled rPython, scLVM, scVEGs (file "scVEGs.R" should be in same directory as "HiVaGe" function), SIEVE
  library(scran)
  library(DropletUtils)
  library(scRNAseq)
  library(magrittr)
  library(Seurat)
  library(scater)
  library(scDblFinder)
  sce_object <- logNormCounts(sce_object)
  valid_flavours <- c("BASiCS", "M3Drop","M3Drop_Basic", "M3Drop_Brennecke", "ROGUE", "ROGUE_n", "Seurat_vst", "Seurat_sct", "Seurat_disp", "scVEGs", "SCHS", "scmap", "SIEVE_Scmap", "SIEVE_Scran", "SIEVE_ROGUE", "SIEVE_M3Drop", "SIEVE_Seurat_vst", "SIEVE_Seurat_disp", "scLVM_counts", "scLVM_log", "scLVM_logvar", "M3Drop_Brennecke_ERCCs", "scLVM_log_ERCCs", "scLVM_logvar_ERCCs", "scLVM_counts_ERCCs", "BASiCS_ERCCs")
  ERCCs_flavours <- c("M3Drop_Brennecke_ERCCs", "scLVM_log_ERCCs", "scLVM_logvar_ERCCs", "scLVM_counts_ERCCs", "BASiCS_ERCCs")
  if (!(flavour %in% valid_flavours)) {
    stop(paste("Invalid flavour. Allowed values are:", paste(valid_flavours, collapse = ", ")))
  } else {
    if ((flavour %in% ERCCs_flavours) && (is.null(ERCCs))){
      stop(paste0("ERRCs not provided for the flavour ", flavour, "."))
    }
  }
  
  if (flavour == "BASiCS") {
    # BASiCS
    # Could use spike-ins
    library(BASiCS)  # Load the BASiCS package if not already loaded
    BASiCS <- newBASiCS_Data(counts(sce_object), Tech = FALSE, SpikeInfo = NULL)
    Chain <- BASiCS_MCMC(Data = BASiCS, N = 1000, Thin = 10, Burn = 500, WithSpikes = FALSE, Regression = TRUE, PrintProgress = FALSE)
    HVGs_BASiCS <- BASiCS_DetectVG(Chain, Task = c("HVG"), PercentileThreshold = round(num_HVGs/length(rownames(sce_object))), VarThreshold = NULL, ProbThreshold = 0.5, EpsilonThreshold = NULL, EFDR = 0.1, Plot = FALSE, MinESS = 100)
    HVGs_BASiCS <- unlist(as.data.frame(HVGs_BASiCS)$GeneName)
    return(HVGs_BASiCS)
    
  } else if (flavour == "M3Drop") {
    # M3Drop 
    # Raw UMI counts
    # Variance-based Feature Selection. Ranks genes by residual dispersion from mean-dispersion power-law relationship.
    # Uses linear regression on log-transformed mean expression and fitted dispersions. Ranks genes by the residuals of this fit, negative = high variance, positive = low variance.
    library(M3Drop)
    HVGs_M3Drop <- sce_object %>%
      counts() %>%
      NBumiConvertData() %>%
      NBumiFitModel() %>%
      NBumiFeatureSelectionHighVar() %>%
      as.data.frame() %>%
      rownames() %>%
      unlist() %>%
      .[1:num_HVGs]
    return(HVGs_M3Drop)
    
  } else if (flavour == "M3Drop_Basic") {
    # M3Drop 
    # Raw UMI counts
    library(M3Drop)
    HVGs_M3Drop <- sce_object %>%
      counts() %>%
      NBumiConvertData() %>%
      NBumiFitBasicModel() %>%
      NBumiFeatureSelectionHighVar() %>%
      as.data.frame() %>%
      rownames() %>%
      unlist() %>%
      .[1:num_HVGs]
    
    return(HVGs_M3Drop)
    
  } else if (flavour == "M3Drop_Brennecke") {
    # Spike-ins is possible
    # Normalized or raw (not log-transformed) expression values, columns = samples, rows = genes.
    library(M3Drop)
    HVGs_M3Drop_Brennecke <- sce_object %>%
      counts() %>%
      BrenneckeGetVariableGenes() %>%
      rownames()
    return(HVGs_M3Drop_Brennecke)
    
  } else if (flavour == "ROGUE") {
    # ROGUE
    # Raw counts
    library(ROGUE)
    count_matrix <- assay(sce_object, "counts") %>%
      as.matrix()
    rownames(count_matrix) <- rownames(sce_object)
    colnames(count_matrix) <- colnames(sce_object)
    ds_ROGUE <- SE_fun(count_matrix, span = 0.5, r = 1, mt.method = "fdr", if.adj = T)
    HVGs_ROGUE <- unlist(ds_ROGUE$Gene)[1:num_HVGs]
    return(HVGs_ROGUE)
    
  } else if (flavour == "ROGUE_n") {
    # ROGUE_n
    # Normalized counts
    library(ROGUE)
    count_matrix <- assay(sce_object, "logcounts") %>%
      as.matrix()
    rownames(count_matrix) <- rownames(sce_object)
    colnames(count_matrix) <- colnames(sce_object)
    ds_ROGUE <- SE_fun(count_matrix, span = 0.5, r = 1, mt.method = "fdr", if.adj = T)
    HVGs_ROGUE_n <- unlist(ds_ROGUE$Gene)[1:num_HVGs]
    return(HVGs_ROGUE_n)
    
  } else if (flavour == "Seurat_vst") {
    # Seurat_vst
    ds_seurat <- sce_object %>%
      logcounts() %>%
      as.data.frame() %>%
      CreateSeuratObject() %>%
      ScaleData() %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = num_HVGs, verbose = TRUE)
    HVGs_seurat_vst <- ds_seurat@assays$RNA@var.features
    return(HVGs_seurat_vst)
    
  } else if (flavour == "Seurat_sct") {
    # Seurat_sct
    ds_seurat <- sce_object %>%
      logcounts() %>%
      as.data.frame() %>%
      CreateSeuratObject() %>%
      ScaleData() %>%
      SCTransform(assay="RNA",variable.features.n = num_HVGs)
    HVGs_seurat_sct <- ds_seurat@assays$SCT@var.features
    return(HVGs_seurat_sct)
    
  } else if (flavour == "Seurat_disp") {
    # Seurat_disp
    ds_seurat <- sce_object %>%
      logcounts() %>%
      as.data.frame() %>%
      CreateSeuratObject() %>%
      ScaleData() %>%
      FindVariableFeatures(selection.method = "disp", nfeatures = num_HVGs, verbose = TRUE)
    HVGs_seurat_disp <- ds_seurat@assays$RNA@var.features
    return(HVGs_seurat_disp)
    
  } else if (flavour == "scVEGs") {
    source("scVEGs.R")
    ds_counts_df <- sce_object %>% 
      counts() %>%
      as.data.frame()
    for (col in 1:ncol(ds_counts_df)) {
      ds_counts_df[, col] <- as.numeric(ds_counts_df[, col])
    }
    pVal <- 0.01
    pFlag <- 1
    species <- 'hs'
    cellSum <- apply(ds_counts_df, 2, sum)
    scaleNum <- mean(cellSum)/cellSum
    scaleFactor <- t(kronecker( matrix(1, 1, dim(ds_counts_df)[1]), scaleNum))
    normData <- scaleFactor * ds_counts_df
    outputName <- 'HVGs_scVEGs'
    HVGs_scVEGs <- scVEGs(normData, pVal, pFlag, species, outputName) %>%
      rownames()
    return(HVGs_scVEGs)
    
  } else if (flavour == "SCHS") {
    library("singleCellHaystack")
    grid_ = 100
    if (ncol(logcounts(sce_object)) < grid_) {
      grid_ = ncol(counts(sce_object))
    }
    set.seed(42)
    HVGs_SCHS <- reducedDim(sce_object, "PCA") %>%
      haystack(expression = logcounts(sce_object), grid.points = grid_) %>%
      show_result_haystack(n = num_HVGs) %>%
      rownames()
    return(HVGs_SCHS)
    
  } else if (flavour == "scmap") {
    library("scmap")
    rowData(sce_object)$feature_symbol <- rownames(sce_object)
    sce_object <- selectFeatures(sce_object, n_features = num_HVGs,
                                 suppress_plot = TRUE) 
    HVGs_scmap <- rowData(sce_object)[rowData(sce_object)$scmap_features == TRUE, ]
    return(rownames(HVGs_scmap))
    
  } else if (flavour == "SIEVE_Scmap") {
    # If your data is a raw counts matrix, method should be chosen from "M3Drop","ROGUE","Scran","Scmap". If your data is a normalized matrix, method should be chosen from "schs"(singleCellHaystack),"Seurat_vst","Seurat_disp","Seurat_sct","ROGUE".
    # SCHS, Seurat_sct doesn't work
    library(SIEVE)
    ds_counts_df <- sce_object %>% 
      counts() %>%
      as.data.frame()
    for (col in 1:ncol(ds_counts_df)) {
      ds_counts_df[, col] <- as.numeric(ds_counts_df[, col])
    }
    ds_counts_seurat <- CreateSeuratObject(counts = ds_counts_df)
    set.seed(42)
    HVGs_SIEVE_Scmap <- fetch_cells(ds_counts_seurat,ratio=0.7,n=50) %>%
      fetch_HVGs_seurat(ds_counts_seurat, ., method="Scmap",n=num_HVGs) %>%
      SIEVE(n=50)
    return(HVGs_SIEVE_Scmap)
    
  } else if (flavour == "SIEVE_Scran") {
    library(SIEVE)
    ds_counts_df <- sce_object %>%
      counts() %>%
      as.data.frame()
    for (col in 1:ncol(ds_counts_df)) {
      ds_counts_df[, col] <- as.numeric(ds_counts_df[, col])
    }
    ds_counts_seurat <- CreateSeuratObject(counts = ds_counts_df)
    set.seed(42)
    HVGs_SIEVE_Scran <- fetch_cells(ds_counts_seurat,ratio=0.7,n=50) %>%
      fetch_HVGs_seurat(ds_counts_seurat, ., method="Scran",n=num_HVGs) %>%
      SIEVE(n=50)
    return (HVGs_SIEVE_Scran)
    
  } else if (flavour == "SIEVE_ROGUE") {
    library(SIEVE)
    ds_counts_df <- sce_object %>%
      logcounts() %>%
      as.data.frame()
    for (col in 1:ncol(ds_counts_df)) {
      ds_counts_df[, col] <- as.numeric(ds_counts_df[, col])
    }
    ds_counts_seurat <- CreateSeuratObject(counts = ds_counts_df)
    set.seed(42)
    HVGs_SIEVE_ROGUE <- fetch_cells(ds_counts_seurat,ratio=0.7,n=50) %>%
      fetch_HVGs_seurat(ds_counts_seurat, ., method="ROGUE",n=num_HVGs) %>%
      SIEVE(n=50)
    return(HVGs_SIEVE_ROGUE)
    
  } else if (flavour == "SIEVE_M3Drop") {
    library(SIEVE)
    ds_counts_df <- sce_object %>% 
      counts() %>%
      as.data.frame()
    for (col in 1:ncol(ds_counts_df)) {
      ds_counts_df[, col] <- as.numeric(ds_counts_df[, col])
    }
    ds_counts_seurat <- CreateSeuratObject(counts = ds_counts_df)
    set.seed(42)
    HVGs_SIEVE_M3Drop <- fetch_cells(ds_counts_seurat,ratio=0.7,n=50) %>%
      fetch_HVGs_seurat(ds_counts_seurat, ., method="M3Drop",n=num_HVGs) %>%
      SIEVE(n=50)
    return(HVGs_SIEVE_M3Drop)
    
  } else if (flavour == "SIEVE_Seurat_vst") {
    library(SIEVE)
    ds_counts_df <- sce_object %>%
      logcounts() %>%
      as.data.frame()
    for (col in 1:ncol(ds_counts_df)) {
      ds_counts_df[, col] <- as.numeric(ds_counts_df[, col])
    }
    ds_counts_seurat <- CreateSeuratObject(counts = ds_counts_df)
    set.seed(42)
    HVGs_SIEVE_Seurat_vst <- fetch_cells(ds_counts_seurat,ratio=0.7,n=50) %>%
      fetch_HVGs_seurat(ds_counts_seurat, ., method="Seurat_vst",n=num_HVGs) %>%
      SIEVE(n=50)
    return(HVGs_SIEVE_Seurat_vst)
    
  } else if (flavour == "SIEVE_Seurat_disp") {
    library(SIEVE)
    ds_counts_df <- sce_object %>%
      logcounts() %>%
      as.data.frame()
    for (col in 1:ncol(ds_counts_df)) {
      ds_counts_df[, col] <- as.numeric(ds_counts_df[, col])
    }
    ds_counts_seurat <- CreateSeuratObject(counts = ds_counts_df)
    set.seed(42)
    HVGs_SIEVE_Seurat_disp <- fetch_cells(ds_counts_seurat,ratio=0.7,n=50) %>%
      fetch_HVGs_seurat(ds_counts_seurat, ., method="Seurat_disp",n=num_HVGs) %>%
      SIEVE(n=50)
    return(HVGs_SIEVE_Seurat_disp)
    
  } else if (flavour == "scLVM_counts") {
    library(scLVM)
    norm_counts <- counts(sce_object) / colData(sce_object)$sizeFactor
    techNoise = fitTechnicalNoise(norm_counts, fit_type = 'counts', use_ERCC = FALSE, plot=FALSE) 
    is_hetLog = getVariableGenes(norm_counts, techNoise$fit, plot=FALSE)
    HVGs_scLVM_counts <- is_hetLog[is_hetLog] %>%
      which() %>%
      names()
    return(HVGs_scLVM_counts)
    
  } else if (flavour == "scLVM_log") {
    library(scLVM)
    # could use spike-ins and size Factors. As threshhold used default
    norm_counts <- counts(sce_object) / colData(sce_object)$sizeFactor
    techNoise = fitTechnicalNoise(norm_counts, fit_type = 'log', use_ERCC = FALSE, plot=FALSE) 
    
    is_hetLog = getVariableGenes(norm_counts, techNoise$fit, plot=FALSE)
    HVGs_scLVM_log <- as.data.frame(is_hetLog[is_hetLog == TRUE]) %>%
      rownames()
    return(HVGs_scLVM_log)
    
  } else if (flavour == "scLVM_logvar") {
    library(scLVM)
    norm_counts <- counts(sce_object) / colData(sce_object)$sizeFactor
    techNoiseLogFit = fitTechnicalNoise(norm_counts, fit_type = 'logvar', use_ERCC = FALSE, plot=FALSE) 
    
    is_hetLog = getVariableGenes(norm_counts, techNoiseLogFit$fit, plot=FALSE)
    HVGs_scLVM_logvar <- as.data.frame(is_hetLog[is_hetLog == TRUE]) %>%
      rownames()
    return(HVGs_scLVM_logvar)
    
  } else if (flavour == "BASiCS_ERCCs") {
    # BASiCS
    # Could use spike-ins
    library(BASiCS)  # Load the BASiCS package if not already loaded
    all_zero_rows <- apply(ERCCs, 1, function(row) all(row == 0))
    ERCCs <- ERCCs[!all_zero_rows, ]
    merged_with_spikes <- rbind(counts(sce_object), ERCCs) %>%
      as.matrix()
    Tech <- merged_with_spikes %>%
      rownames() %>%
      grepl("ERCC",.)
    SpikeInfo <- rowSums(ERCCs[, -1]) %>%
      data.frame(Genes = rownames(ERCCs), ERCCs_Sum = .)
    BASiCS <- newBASiCS_Data(merged_with_spikes, Tech = Tech, SpikeInfo = SpikeInfo)
    Chain <- BASiCS_MCMC(Data = BASiCS, N = 1000, Thin = 10, Burn = 500, WithSpikes = TRUE, Regression = TRUE, PrintProgress = FALSE)
    HVGs_BASiCS <- BASiCS_DetectVG(Chain, Task = c("HVG"), PercentileThreshold = round(num_HVGs/length(rownames(sce_object))), VarThreshold = NULL, ProbThreshold = 0.5, EpsilonThreshold = NULL, EFDR = 0.1, Plot = FALSE, MinESS = 100)
    HVGs_BASiCS <- unlist(as.data.frame(HVGs_BASiCS)$GeneName)
    
    return(BASiCS)
    
  }
  else if (flavour == "scLVM_counts_ERCCs") {
    library(scLVM)
    norm_counts <- counts(sce_object) / colData(sce_object)$sizeFactor
    techNoise = fitTechnicalNoise(norm_counts, nCountsERCC = ERCCs, fit_type = 'counts', use_ERCC = TRUE, plot=FALSE) 
    is_hetLog = getVariableGenes(norm_counts, techNoise$fit, plot=FALSE)
    HVGs_scLVM_counts <- is_hetLog[is_hetLog] %>%
      which() %>%
      names()
    return(HVGs_scLVM_counts)
    
  } else if (flavour == "scLVM_log_ERCCs") {
    library(scLVM)
    norm_counts <- counts(sce_object) / colData(sce_object)$sizeFactor
    techNoise = fitTechnicalNoise(norm_counts, nCountsERCC = ERCCs, fit_type = 'log', use_ERCC = TRUE, plot=FALSE) 
    is_hetLog = getVariableGenes(norm_counts, techNoise$fit, plot=FALSE)
    HVGs_scLVM_log <- is_hetLog[is_hetLog] %>%
      which() %>%
      names()
    return(HVGs_scLVM_log)
    
  } else if (flavour == "scLVM_logvar_ERCCs") {
    library(scLVM)
    norm_counts <- counts(sce_object) / colData(sce_object)$sizeFactor
    techNoise = fitTechnicalNoise(norm_counts, nCountsERCC = ERCCs, fit_type = 'logvar', use_ERCC = TRUE, plot=FALSE) 
    is_hetLog = getVariableGenes(norm_counts, techNoise$fit, plot=FALSE)
    HVGs_scLVM_logvar <- is_hetLog[is_hetLog] %>%
      which() %>%
      names()
    return(HVGs_scLVM_logvar)
    
  } else if (flavour == "M3Drop_Brennecke_ERCCs") {
    # Spike-ins is possible
    # Normalized or raw (not log-transformed) expression values, columns = samples, rows = genes.
    library(M3Drop)
    HVGs_M3Drop_Brennecke <- sce_object %>%
      counts() %>%
      BrenneckeGetVariableGenes(., ERCCs) %>%
      rownames()
    return(HVGs_M3Drop_Brennecke)
    
  } 
}


### Main function
# sce_object - SingleCellExperiment variable
# labels - vector with types of the cells or reference-annotated cell types
# batch, vector containing batch of each cell
# num_HVGs = 400, num of HVGs to get
# assay.type = c("counts" or "logcounts"), sce_object's assay to use for some metrics
HiVaGeMetrics <- function(sce_object, labels = NULL, batch = NULL, num_HVGs = 400, assay.type = "counts") {
  # Methods of getting HVGs
  R_flavours = c("M3Drop","M3Drop_Basic", "M3Drop_Brennecke", "ROGUE", "ROGUE_n", "Seurat_vst", "Seurat_sct", "Seurat_disp", "scVEGs", "SCHS", "scmap", "SIEVE_Scmap", "SIEVE_Scran", "SIEVE_ROGUE", "SIEVE_M3Drop", "SIEVE_Seurat_vst", "SIEVE_Seurat_disp", "scLVM_counts", "scLVM_log", "scLVM_logvar", "BASiCS")
  Py_flavours = c('scanpy_seurat', 'scanpy_cell_ranger', 'scanpy_Pearson',  "scanpy_seurat_v3", 'Triku')
  R_ERCCs_flavours <- c("M3Drop_Brennecke_ERCCs", "scLVM_log_ERCCs", "scLVM_logvar_ERCCs", "scLVM_counts_ERCCs", "BASiCS_ERCCs")
  # Full range
  valid_flavours <- c(R_flavours, Py_flavours)
  # Flavours grouped by language
  flavour_flavours = list("R_flavours" = R_flavours, "Py_flavours" = Py_flavours)
  # If ERCCs, add more methods
  if (altExpNames(sce_object)[1] == "ERCC") {
    valid_flavours = c(valid_flavours, R_ERCCs_flavours)
    flavour_flavours[["R_ERCCs_flavours"]] = R_ERCCs_flavours
  }
  
  # Prepare results dataframe
  metrics = c("Overlap with highly expressed genes", "Overlap with lowly expressed genes",
              "Pearson's correlation between mean and variance","Calinski-Harabasz index",
              "Davies-Bouldin index","Average silhouette width", "Adjusted Rand index",
              "Average ROGUE score", "Purity score", "Runtime", "hvgs")
  res_df = data.frame(matrix(ncol = length(metrics))) 
  colnames(res_df) = metrics
  # Prepare list with other metrics
  other_metrics = c("rogue_score_boxplot", "hvgs_correlation")
  other_res_list = list()
  for (i in valid_flavours) {
    other_res_list[[i]] = NA
  }

  # Counter for metrics
  k = 1

  # Name of the folder to save data for
  hvgs_folder = "NAME"
  dir.create(paste0(hvgs_folder,"_HVGs_csv"), showWarnings = FALSE)

  # Loop through flavours
  for (group in names(flavour_flavours)) {
    
    if (group == "R_ERCCs_flavours"){
      # Look for ERCCs in SingleCellExperiment
      replaceCat("Preparing ERCCs variable.. \n")
      ERCCs <- sce_object |>
        altExp() |>
        counts()
    } else if (group == "Py_flavours") {
      replaceCat("Preparing Python script and variables.. \n")
      # Create dataframe for python methods
      py_df = as.data.frame(counts(sce_object))
      py_df$ProbeIDs = row.names(py_df)
      py_df <- py_df |>
        dplyr::select(ProbeIDs, everything())
      # Convert to pandas dataframe
      py$ds = r_to_py(py_df)
      # Source HVGs script (supply path to python methods script)
      #source_python("HiVaGePY.py")
    }
    
    for (i in flavour_flavours[[group]]) {
      # Prepare results variable for this specific method
      res = list()
      for (j in metrics){
        res[[j]] = NA
      }
      
      res_other = list()
      for (j in other_metrics){
        res_other[[j]] = NA
      }
      
      # Get HVGs
      replaceCat(paste("Starting ", i, ".\n", sep = ""))
      replaceCat(paste("Getting HVGs (", i, ")..", sep = ""))

      # Calculate runtime for the method
      if (i %in% R_flavours) {
        tic()
        hvgs = HiVaGe(sce_object, i, batching = batch, num_HVGs = 400)
        time_temp = toc()
        res[["Runtime"]] = as.numeric(time_temp$toc-time_temp$tic)
      } else if (i %in% Py_flavours) {
        tic()
        hvgs = HiVaGePY(py$ds, i, HVGs_num = 400)$ProbeIDs
        time_temp = toc()
        res[["Runtime"]] = as.numeric(time_temp$toc-time_temp$tic)
      } else if (i %in% R_ERCCs_flavours) {
        tic()
        hvgs = HiVaGe(sce_object, i, batching = batch, num_HVGs = 400, ERCCs)
        time_temp = toc()
        res[["Runtime"]] = as.numeric(time_temp$toc-time_temp$tic)
      } else {
        replaceCat("Invalid flavor \"", i, "\".\n", sep = "")
      }
      
      # Checking if at least 1 HVG found
      if (length(hvgs) > 1){
        # Saving list of found HVGs
        res[["hvgs"]] = paste(hvgs, collapse = ", ")
        
        replaceCat("##Clustering.\n")
        set.seed(42)
        ## Clustering based on HVGs
        # Get data from sce_objectExperiment for the HVGs
        filtered <- sce_object[hvgs, ]
        # GLMPCA should reduce dimensions count to amount less then there are HVGs
        L_ = 10
        if (length(hvgs) < 10){
          L_ = round(length(hvgs)/2)
          replaceCat(paste("Low HVG count (< 10), reducing dimensions to (", length(hvgs), "/2) = ", L_, " dimensions!\n", sep = ""))
        }
        replaceCat("Running GLMPCA for dimensionality reduction..")
        filtered <- GLMPCA(filtered, L=L_, minibatch="stochastic")
        replaceCat("Building shared nearest-neighbor graph (SNNGraph)..") 
        g <- buildSNNGraph(filtered, k=10, use.dimred = 'GLMPCA')
        replaceCat("Performing Louvain clustering..")
        clust <- cluster_louvain(g)
        # Add cluster data as a factor to "filtered" variable
        filtered$Louvain <- factor(membership(clust))
        
        replaceCat("##Moving onto metrics.\n")
        ### Clustering metrics
        # Get relevant data
        if (assay.type == "counts"){
          sce_data = t(counts(filtered))
        } else if (assay.type == "logcounts"){
          sce_data = t(logcounts(filtered))
        } else {
          stop("\nInnappropriate assay type.\n")
        }
        # Make sure clustering vector is proper data type
        clustering = as.integer(factor(filtered$Louvain))
        
        replaceCat("Calculating CH and DB indexes..")
        res[["Calinski-Harabasz index"]] = index.G1(sce_data, clustering)
        res[["Davies-Bouldin index"]] = index.DB(sce_data, clustering)$DB
        
        replaceCat("Calculating average silhouette width.. (this one takes a bit)")
        # First get distance matrix
        sce_dist_matrix = dist(sce_data)
        # Calculate silhouettes
        si_asw = silhouette(clustering, sce_dist_matrix)
        # Get average silhouette width
        res[["Average silhouette width"]] = summary(si_asw)$avg.width
        
        if (!is.null(labels)) {
          replaceCat("Calculating Purity score.. (this one takes a bit)\n")
          res[["Purity score"]] = getAvrPurity(hvgs, sce_object, labels, tSNE_count = 3)
          
          # Make sure labels vector is proper data type for next metric
          labels = as.integer(factor(labels))
          replaceCat("Calculating Adjusted Rand index..")
          res[["Adjusted Rand index"]] = adj.rand.index(labels, clustering)
        } else {
          replaceCat("\nSkipping Purity score and Adjusted Rand Index since labels vector not provided.\n")
        }
        
        if (!is.null(batch)) {
          # ROGUE score is prone to breaking, so it is turned off
          #replaceCat("Calculating ROGUE score..\n")
          #temp_rogue = getROGUEScore(hvgs, filtered, clustering = filtered$Louvain,
          #                           sampling = batch, platform = "full-length",
          #                           assay.type = assay.type)
          
          #res[["Average ROGUE score"]] = mean(apply(temp_rogue, 1, mean, na.rm=TRUE), na.rm=TRUE)
          #res_other[["rogue_score_boxplot"]] = rogue.boxplot(as.data.frame(t(temp_rogue)))
          
          replaceCat("Calculating dependency with mean expression..")
          depend_list = getOverlapWithHighLowExpressed(hvgs, sce_object, batching = batch)
          res[["Overlap with highly expressed genes"]] = depend_list$highly_overlap
          res[["Overlap with lowly expressed genes"]] = depend_list$lowly_overlap
          res[["Pearson's correlation between mean and variance"]] = depend_list$correlation_meanVariance
        } else {
          replaceCat("\nSkipping ROGUE score and Dependency with mean expression since batch vector not provided.\n")
        }
        
        replaceCat("Calculating correlation between HVGs (this one takes a bit)..")
        res_other[["hvgs_correlation"]] = correlatedHVGs(hvgs, filtered)

        # Updating results dataframe
        res_df = rbind(res_df, res)
        k = k + 1
        row.names(res_df)[k] = i
        # Same for other metrics
        other_res_list[[i]] = res_other

        # Save all metrics together
        full_res = list(metrics_df = res_df, other = other_res_list)
        # to a temporary file in case a crash happens
        saveRDS(full_res, file = file.path(paste0(hvgs_folder,"_HVGs_csv"), paste0(i, "_", hvgs_folder, "_backup.rds")))
        replaceCat(paste("##", i, " done.\n", sep = ""))
        
      } else {
        replaceCat(paste("##No HVGs found with ", i, "! Skipping.\n", sep = ""))
      }
    }
  }

  # Save all metrics together properly
  full_res = list(metrics_df = res_df[-1,], other = other_res_list)
  # Return the list with a dataframe and a list with all metrics
  return(full_res)
}
