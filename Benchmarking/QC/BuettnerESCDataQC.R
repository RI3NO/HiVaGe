library(scRNAseq)
library(scater)
library(edgeR)
library(DropletUtils)
library(patchwork)

library(ensembldb)
library(EnsDb.Mmusculus.v75)

library(AnnotationHub)

Buettner_ds = BuettnerESCData()

bcrank <- barcodeRanks(counts(Buettner_ds))
# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2)

abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

View(as.data.frame(bcrank[order(bcrank$total),][1:50,]))

# just look at those last 9 cells, they're miserable
Buettner_ds = Buettner_ds[,rownames(bcrank[order(bcrank$total),][-(1:9),])]


# remove low ave.counts
allzero <- rowMeans(counts(Buettner_ds) == 0) == 1
# Create a table showing the number of cells with all-zero counts
table(allzero)
# Remove cells with all-zero counts for all genes from 'sce'
Buettner_ds <- Buettner_ds[which(!allzero),]

# logcounts for PCA
Buettner_ds <- logNormCounts(Buettner_ds)
# Dimensionality reduction with PCA
# Run PCA on 'sce'
Buettner_ds <- runPCA(Buettner_ds)
# Compute doublet density using the reduced dimensions from PCA
dbl.dens <- computeDoubletDensity(Buettner_ds, d = ncol(reducedDim(Buettner_ds)))
# Add doublet scores to the 'sce' object
Buettner_ds$DoubletScore <- dbl.dens
# Remove doublets from 'sce' based on the identified doublet calls
# Identify doublets using doubletThresholding function
dbl.calls <- doubletThresholding(data.frame(score = dbl.dens),
                                 method = "griffiths", returnType = "call")
# table(dbl.calls == "singlet")
# Remove doublets from 'sce' based on the identified doublet calls
Buettner_ds <- Buettner_ds[, dbl.calls == "singlet"]
# Perform PCA again
Buettner_ds <- runPCA(Buettner_ds)

# not required as per how the knee plot looks + e.out function themselves tells this dataset has already been worked on
#e.out <- emptyDrops(counts(Buettner_ds))
#Buettner_ds <- Buettner_ds[,which(e.out$FDR <= 0.001)]

# rowData(Buettner_ds)$is_spike = grepl("ERCC", rownames(Buettner_ds))
# altExp(Buettner_ds, "ERCC") = Buettner_ds[rowData(Buettner_ds)$is_spike]
# Buettner_ds = Buettner_ds[!rowData(Buettner_ds)$is_spike]

#loading old annotation dataset
library(ensembldb)
library(EnsDb.Mmusculus.v75)

ens.mm.v75 = EnsDb.Mmusculus.v75
chr.loc <- mapIds(ens.mm.v75, keys=rownames(Buettner_ds),
                  keytype="GENEID", column="SEQNAME")
is.mito <- which(chr.loc=="MT")

# which.mito == which(ensambl=="MT")
# threshold_feature can be for example "subsets_Mito_percent" or "altexps_ERCC_percent" or any other QC metric
getDistributionPlots <- function(singleCell, which.mito, threshold_feature = "sum"){
  singleCell <- addPerCellQC(singleCell)
  if (!missing(which.mito)) {
    singleCell <- addPerCellQC(singleCell, subsets=list(Mito=which.mito))
  }
  
  df = perCellQCMetrics(singleCell)
  pos_batches = names(colData(singleCell))
  pos_batches = pos_batches[!pos_batches %in% colnames(df)]
  
  p_1 = list()
  p_2 = list()
  
  for (i in pos_batches) {
    batch_ = singleCell@colData@listData[[i]]
    
    #batch.reasons <- perCellQCFilters(df, batch=batch_,
    #                                  sub.fields=c("altexps_ERCC_sum"))
    #colSums(as.matrix(batch.reasons))
    
    discard.ercc <- isOutlier(singleCell[[threshold_feature]],
                              type="higher", batch=batch_)
    p = plotColData(singleCell, x=i, y=threshold_feature,
                    colour_by=I(discard.ercc))
    p_1[[i]] = p
    
    ercc.thresholds <- attr(discard.ercc, "thresholds")["higher",]
    ercc.thresholds
    names(ercc.thresholds)[isOutlier(ercc.thresholds, type="higher")]
    
    discard.ercc2 <- isOutlier(singleCell[[threshold_feature]],
                               type="higher", batch=batch_,
                               subset=!(batch_ %in% names(ercc.thresholds)[isOutlier(ercc.thresholds, type="higher")]))
    
    p_ = plotColData(singleCell, x=i, y=threshold_feature,
                     colour_by=I(discard.ercc2))
    p_2[[i]] = p_
  }
  
  discard.ercc <- isOutlier(singleCell[[threshold_feature]],
                            type="higher")
  p = plotColData(singleCell, y=threshold_feature,
                  colour_by=I(discard.ercc))
  p_3 = p
  
  return(list(batch_based = p_1, good_batch_based = p_2, general = p_3))
}

plots_ = getDistributionPlots(Buettner_ds, is.mito, "subsets_Mito_percent")
plots_$good_batch_based$phase
# cell phase freaking messing stuff up!
# we only taking one cell phase
Buettner_ds = Buettner_ds[,colData(Buettner_ds)$phase == "G1"]

stats <- perCellQCMetrics(Buettner_ds, subsets=list(Mito=is.mito))
high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
Buettner_ds = Buettner_ds[,!high.mito]

df <- perCellQCMetrics(Buettner_ds)
reasons <- perCellQCFilters(df)
Buettner_ds <- Buettner_ds[,!reasons$discard]

# Load reference data for SingleR
ref <- ImmGenData(ensembl = TRUE)

# Perform SingleR prediction on 'Buettner_ds' using the reference data
pred <- SingleR(test = Buettner_ds, ref = ref, labels = ref$label.fine, assay.type.test = "logcounts")

# Assign the predicted cell types to 'Buettner_ds$CellTypes'
Buettner_ds$CellTypes <- factor(pred$pruned.labels)

Buettner_ds$CellTypes
