# Loading ds, libraries 
library(scRNAseq)
library(DESeq2)
sce <- RichardTCellData()


ERCCs <- sce_subset %>%
  altExp() %>%
  counts()

selected_genes <- rownames(sce) %in% rownames(sce)[1:200]
sce_subset <- sce[selected_genes, ]

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
    techNoise = fitTechnicalNoise(norm_counts, nCountsERCC = ERCCs, fit_type = 'counts', use_ERCC = TRUE, plot=FALSE) 
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
    techNoise = fitTechnicalNoise(norm_counts, nCountsERCC = CountsERCC, fit_type = 'log', use_ERCC = TRUE, plot=FALSE) 
    is_hetLog = getVariableGenes(norm_counts, techNoise$fit, plot=FALSE)
    HVGs_scLVM_log <- is_hetLog[is_hetLog] %>%
      which() %>%
      names()
    return(HVGs_scLVM_log)
    
  } else if (flavour == "scLVM_logvar_ERCCs") {
    library(scLVM)
    norm_counts <- counts(sce_object) / colData(sce_object)$sizeFactor
    techNoise = fitTechnicalNoise(norm_counts, nCountsERCC = CountsERCC, fit_type = 'logvar', use_ERCC = TRUE, plot=FALSE) 
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


# Example
HVGs_1 <- HiVaGe(sce_subset, "sce", 0.1, ERCCs)
valid_flavours <- c("BASiCS", "M3Drop", "M3Drop_Basic", "M3Drop_Brennecke", "ROGUE", "ROGUE_n", "Seurat_vst", "Seurat_sct", "Seurat_disp", "scVEGs", "SCHS", "scmap", "SIEVE_Scmap", "SIEVE_Scran", "SIEVE_ROGUE", "SIEVE_M3Drop", "SIEVE_Seurat_vst", "SIEVE_Seurat_disp", "scLVM_counts", "scLVM_log", "scLVM_logvar", "M3Drop_Brennecke_ERCCs", "scLVM_log_ERCCs", "scLVM_logvar_ERCCs", "scLVM_counts_ERCCs", "BASiCS_ERCCs")


