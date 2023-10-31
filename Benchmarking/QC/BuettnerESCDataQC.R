library(scRNAseq)
library(scater)
library(edgeR)
library(DropletUtils)
library(patchwork)

library(ensembldb)
library(EnsDb.Mmusculus.v75)

library(AnnotationHub)

# Loading the dataset from scRNAseq library
Buettner_ds = BuettnerESCData()

## Knee plot preparation
bcrank <- barcodeRanks(counts(Buettner_ds))
# Only showing unique points for plotting speed
uniq <- !duplicated(bcrank$rank)
# Knee plot itself
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)
legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

# Look at the last 50 cells based on barcode ranks
View(as.data.frame(bcrank[order(bcrank$total),][1:50,]))
# Just look at those last 9 cells, they're miserable
# Remove them
Buettner_ds = Buettner_ds[,rownames(bcrank[order(bcrank$total),][-(1:9),])]

# Removing empty drops is not required as per how the knee plot looks
# In addition, e.out function itself tells this dataset has already been filtered for empty droplets

# Remove low ave.counts
allzero <- rowMeans(counts(Buettner_ds) == 0) == 1
# Create a table showing the number of cells with all-zero counts
table(allzero)
# Remove cells with all-zero counts for all genes from 'sce'
Buettner_ds <- Buettner_ds[which(!allzero),]

# Logcounts for PCA
Buettner_ds <- logNormCounts(Buettner_ds)
# Dimensionality reduction with PCA
# Run PCA
Buettner_ds <- runPCA(Buettner_ds)
# Compute doublet density using the reduced dimensions from PCA
dbl.dens <- computeDoubletDensity(Buettner_ds, d = ncol(reducedDim(Buettner_ds)))
# Add doublet scores to the SingleCellExperiment object
Buettner_ds$DoubletScore <- dbl.dens
# Remove doublets from SingleCellExperiment based on the identified doublet calls
# Identify doublets using doubletThresholding function
dbl.calls <- doubletThresholding(data.frame(score = dbl.dens),
                                 method = "griffiths", returnType = "call")
# Remove doublets from SingleCellExperiment based on the identified doublet calls
Buettner_ds <- Buettner_ds[, dbl.calls == "singlet"]
# Perform PCA again
Buettner_ds <- runPCA(Buettner_ds)

# Loading the annotation dataset
library(EnsDb.Mmusculus.v75)
ens.mm.v75 = EnsDb.Mmusculus.v75
# Get a vector containing gene loci information
chr.loc <- mapIds(ens.mm.v75, keys=rownames(Buettner_ds),
                  keytype="GENEID", column="SEQNAME")
# Vector containing information about if gene is mitochondrial
is.mito <- which(chr.loc=="MT")

# which.mito == which(ensambl=="MT")
# threshold_feature can be for example "subsets_Mito_percent" or "altexps_ERCC_percent" or any other QC metric
getDistributionPlots <- function(singleCell, which.mito, threshold_feature = "sum"){
  singleCell <- addPerCellQC(singleCell)
  # Get QC mito metric if it's missing 
  if (!missing(which.mito)) {
    singleCell <- addPerCellQC(singleCell, subsets=list(Mito=which.mito))
  }

  # Get all QC metrics
  df = perCellQCMetrics(singleCell)
  # Get possible batches from metadata
  pos_batches = names(colData(singleCell))
  # Remove the added QC metrics from them
  pos_batches = pos_batches[!pos_batches %in% colnames(df)]

  # Prepare result lists
  p_1 = list()
  p_2 = list()

  # Loop through metadata columns
  for (i in pos_batches) {
    # Temporary batch value
    batch_ = singleCell@colData@listData[[i]]

    # Get outliers based on this temp batch
    discard.ercc <- isOutlier(singleCell[[threshold_feature]],
                              type="higher", batch=batch_)
    # Make the violin plot
    p = plotColData(singleCell, x=i, y=threshold_feature,
                    colour_by=I(discard.ercc))
    # Save the plot
    p_1[[i]] = p

    # Look for batches that are outliers themselves
    ercc.thresholds <- attr(discard.ercc, "thresholds")["higher",]
    # Remove them from the outlier detection
    discard.ercc2 <- isOutlier(singleCell[[threshold_feature]],
                               type="higher", batch=batch_,
                               subset=!(batch_ %in% names(ercc.thresholds)[isOutlier(ercc.thresholds, type="higher")]))
    # New violin plot that ignores outlier batches
    p_ = plotColData(singleCell, x=i, y=threshold_feature,
                     colour_by=I(discard.ercc2))
    # Save it
    p_2[[i]] = p_
  }

  # Discard outliers without separation of cells between batches
  discard.ercc <- isOutlier(singleCell[[threshold_feature]],
                            type="higher")
  # Violin plot with no batches
  p = plotColData(singleCell, y=threshold_feature,
                  colour_by=I(discard.ercc))
  # Save it
  p_3 = p

  # Return all plots
  return(list(batch_based = p_1, good_batch_based = p_2, general = p_3))
}

# Get different violin plots to check for which variables need to be used for as different batches in outlier detection
plots_ = getDistributionPlots(Buettner_ds, is.mito, "subsets_Mito_percent")
plots_$good_batch_based$phase
# Cell phase messing stuff up!
# We only taking one cell phase
Buettner_ds = Buettner_ds[,colData(Buettner_ds)$phase == "G1"]

# Get QC metrics based on mito genes
stats <- perCellQCMetrics(Buettner_ds, subsets=list(Mito=is.mito))
# Cells with high mito are outliers
high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
# Remove said cells
Buettner_ds = Buettner_ds[,!high.mito]

# Get QC the other metrics
df <- perCellQCMetrics(Buettner_ds)
reasons <- perCellQCFilters(df)
# Remove outliars
Buettner_ds <- Buettner_ds[,!reasons$discard]

# Load reference data for SingleR
ref <- ImmGenData(ensembl = TRUE)

# Perform SingleR prediction on 'Buettner_ds' using the reference data
pred <- SingleR(test = Buettner_ds, ref = ref, labels = ref$label.fine, assay.type.test = "logcounts")

# Assign the predicted cell types to 'Buettner_ds$CellTypes'
Buettner_ds$CellTypes <- factor(pred$pruned.labels)
