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
source_python("R/flowers_and_clouds/HiVaGePY.py")

library(magrittr)
library(Seurat)
source("scVEGs.R")
library("singleCellHaystack")
library("scmap")
library(SIEVE)


BiocManager::install("RPushbullet")
library(RPushbullet)
pbSetup()





test_percentile = function(sce_object, flavour, labels, batching = 1, step = 0.1, start = 0.05, repeats = 3) {
  R_flavours = c("scmap", "M3Drop", "M3Drop_Basic", "M3Drop_Brennecke", "ROGUE", "ROGUE_n", "Seurat_vst", "Seurat_sct", "Seurat_disp", "scVEGs", "SCHS", "SIEVE_Scmap", "SIEVE_Scran", "SIEVE_ROGUE", "SIEVE_M3Drop", "SIEVE_Seurat_vst", "SIEVE_Seurat_disp") #"BASiCS") #"scLVM_log", "scLVM_logvar")
  Py_flavours = c('scanpy_seurat', 'scanpy_cell_ranger', 'scanpy_Pearson', 'Triku') #"scanpy_seurat_v3"  #
  
  metrics = c("Calinski-Harabasz index", "Davies-Bouldin index", "Adjusted Rand index", "percentile", "num_HVGs")
  res_df = data.frame(matrix(ncol = length(metrics)))
  colnames(res_df) = metrics
  k = 1
  
  
  #replaceCat("Preparing Python script and variables..")
  # Create dataframe for python methods
  # py_df = as.data.frame(counts(sce))
  # py_df$ProbeIDs = row.names(py_df)
  # py_df <- py_df |>
  #   dplyr::select(ProbeIDs, everything())
  # # Convert to pandas dataframe
  # py$ds = r_to_py(py_df)
  # # Source HVGs script (supply path to python methods script)
  # #source_python("HiVaGePY.py")
  
  
  replaceCat(paste("### Start percentile: ", start, "\n Step: 0.1", "\n Repeats:", repeats, ".\n", sep = ""))
  
  for (i in 0:repeats){
    perc = step*i+start
    replaceCat(paste("##Starting with percentile = ", perc, ". Repeat #",i," out of ",repeats,".\n", sep = ""))
    
    if (flavour %in% R_flavours) {
      hvgs = HiVaGe(sce_object, flavour, batching = batching, percentile = perc)
    } else if (flavour %in% Py_flavours) {
      hvgs = HiVaGePY(py$ds, flavour, perc)$ProbeIDs
    } else {
      replaceCat("Invalid flavor \"", i, "\".\n", sep = "")
    }
    
    replaceCat("##Clustering.\n")
    set.seed(42)
    ## Clustering based on HVGs
    # Get data from sce_objectExperiment for the HVGs
    if (length(hvgs) > 0) {
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
      # filtered <- runTSNE(filtered, dimred="GLMPCA")
      # plotTSNE(filtered, colour_by="Louvain")
      
      replaceCat("##Moving onto metrics.\n")
      # Prepare results variable
      res = list()
      for (j in metrics){
        res[[j]] = NA
      }
      
      ### Clustering metrics
      # Get relevant data
      sce_data = t(counts(sce_object))
      
      # Make sure clustering vector is proper data type
      clustering = as.integer(factor(filtered$Louvain))
      
      replaceCat("Calculating CH and DB indexes..")
      res[["Calinski-Harabasz index"]] = index.G1(sce_data, clustering)
      res[["Davies-Bouldin index"]] = index.DB(sce_data, clustering)$DB
      if (!missing(labels)) {
        # Make sure labels vector is proper data type for next metric
        labels = as.integer(factor(labels))
        replaceCat("Calculating Adjusted Rand index..")
        res[["Adjusted Rand index"]] = adj.rand.index(labels, clustering)
      } else {
        replaceCat("\nSkipping Adjusted Rand Index since labels vector not provided.\n")
      }
    } else {
      replaceCat("##No HVGs found!\n")
    }
    
    res[["percentile"]] = as.numeric(perc)
    res[["num_HVGs"]] = (length(rownames(sce_object)) * perc)
    res_df = rbind(res_df, res)
    
    
    k = k + 1
    row.names(res_df)[k] = perc
  }
  
  return(res_df[-1,])
}





Aztekin_ds_scmap = try(test_percentile(Aztekin_ds, "scmap", labels = Aztekin_ds$cluster, batching = Aztekin_ds$new_batch, step = 0.1, start = 0.05, repeats = 7))
if(inherits(Aztekin_ds_scmap, "try-error")){
  pbPost("note", "R", "Error encountered. scmap")
} else {
  Aztekin_ds_SCHS = try(test_percentile(Aztekin_ds, "SCHS", labels = Aztekin_ds$cluster, batching = Aztekin_ds$new_batch, step = 0.1, start = 0.05, repeats = 7))
  if(inherits(Aztekin_ds_SCHS, "try-error")){
    pbPost("note", "R", "Error encountered. SCHS")
  } else {
    Aztekin_ds_Seurat_disp = try(test_percentile(Aztekin_ds, "Seurat_disp", labels = Aztekin_ds$cluster, batching = Aztekin_ds$new_batch, step = 0.1, start = 0.05, repeats = 7))
    if(inherits(Aztekin_ds_Seurat_disp, "try-error")){
      pbPost("note", "R", "Error encountered. Seurat_disp")
    } else {
      pbPost("note", "R", "0.05-0.75 done!")
      Aztekin_ds_scmap_fine = try(test_percentile(Aztekin_ds, "scmap", labels = Aztekin_ds$cluster, batching = Aztekin_ds$new_batch, step = 0.01, start = 0.005, repeats = 7))
      if(inherits(Aztekin_ds_scmap_fine, "try-error")){
        pbPost("note", "R", "Error encountered. scmap_fine")
      } else {
        Aztekin_ds_SCHS_fine = try(test_percentile(Aztekin_ds, "SCHS", labels = Aztekin_ds$cluster, batching = Aztekin_ds$new_batch, step = 0.01, start = 0.005, repeats = 7))
        if(inherits(Aztekin_ds_SCHS_fine, "try-error")){
          pbPost("note", "R", "Error encountered. SCHS_fine")
        } else {
          Aztekin_ds_Seurat_disp_fine = try(test_percentile(Aztekin_ds, "Seurat_disp", labels = Aztekin_ds$cluster, batching = Aztekin_ds$new_batch, step = 0.01, start = 0.005, repeats = 7))
          if(inherits(Aztekin_ds_Seurat_disp_fine, "try-error")){
            pbPost("note", "R", "Error encountered. Seurat_disp_fine")
          } else {
            pbPost("note", "R", "All done!")
          }
        }
      }
    }
  }
}

library(ggplot2)
colors <- c("scmap" = "blue", "SCHS" = "red", "seurat_disp" = "orange")
ggplot(mapping = aes(x=percentile, y=`Adjusted Rand index`)) + 
  geom_line(data=rbind(Aztekin_ds_scmap, Aztekin_ds_scmap_fine), aes(color='scmap')) +
  geom_line(data=rbind(Aztekin_ds_SCHS, Aztekin_ds_SCHS_fine), aes(color='SCHS')) +
  geom_line(data=rbind(Aztekin_ds_suertat_disp, Aztekin_ds_suertat_disp_fine), aes(color='seurat_disp')) +
  scale_colour_manual(values = colors) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = round(seq(0, 1, by = 0.05),3)) +
  theme_bw()


ggplot(mapping = aes(x=num_HVGs, y=`Adjusted Rand index`)) + 
  geom_line(data=rbind(Aztekin_ds_scmap, Aztekin_ds_scmap_fine), aes(color='scmap')) +
  geom_line(data=rbind(Aztekin_ds_SCHS, Aztekin_ds_SCHS_fine), aes(color='SCHS')) +
  geom_line(data=rbind(Aztekin_ds_suertat_disp, Aztekin_ds_suertat_disp_fine), aes(color='seurat_disp')) +
  scale_colour_manual(values = colors) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 350)) +
  theme_bw()


ggplot(mapping = aes(x=percentile, y=`Calinski-Harabasz index`)) + 
  geom_line(data=rbind(Aztekin_ds_scmap, Aztekin_ds_scmap_fine), aes(color='scmap')) +
  geom_line(data=rbind(Aztekin_ds_SCHS, Aztekin_ds_SCHS_fine), aes(color='SCHS')) +
  geom_line(data=rbind(Aztekin_ds_suertat_disp, Aztekin_ds_suertat_disp_fine), aes(color='seurat_disp')) +
  scale_colour_manual(values = colors) +
  #scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = round(seq(0, 1, by = 0.05),3)) +
  theme_bw()


ggplot(mapping = aes(x=percentile, y=`Davies-Bouldin index`)) + 
  geom_line(data=rbind(Aztekin_ds_scmap, Aztekin_ds_scmap_fine), aes(color='scmap')) +
  geom_line(data=rbind(Aztekin_ds_SCHS, Aztekin_ds_SCHS_fine), aes(color='SCHS')) +
  geom_line(data=rbind(Aztekin_ds_suertat_disp, Aztekin_ds_suertat_disp_fine), aes(color='seurat_disp')) +
  scale_colour_manual(values = colors) +
  #scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = round(seq(0, 1, by = 0.05),3)) +
  theme_bw()