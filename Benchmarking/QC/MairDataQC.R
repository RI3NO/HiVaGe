library(scRNAseq)
library(EnsDb.Hsapiens.v86)
Mair_ds <- MairPBMCData()
View(Mair_ds)

rownames(Mair_ds) <- rowData(Mair_ds)$Symbol
new_rownames = gsub("_", "-", rownames(Mair_ds))
rownames(Mair_ds) = make.names(new_rownames, unique=TRUE)

# Making a plot that shows that there's clear distinction in cells that have big UMIs and low UMIs. Inflection is 100 UMIs so we could use default threshold. Big difference in genes expression we consider normal, as we have PBMC cells, most of which are specific to some function and don't express big amount of genes.
bcrank <- barcodeRanks(counts(Mair_ds))

# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2)

abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

# Deleting empty droplets using false discovery rate threshhold
e.out <- emptyDrops(counts(Mair_ds))
summary(e.out$FDR <= 0.001)
table(colSums(counts(Mair_ds))>100, e.out$FDR<=0.001, useNA = "ifany")
# By default all droplets with fewer than 100 UMIs are considered empty.
Mair_ds <- Mair_ds[, which(e.out$FDR <= 0.001)]

# Delete cells that have 0 genes expressed and genes that are expressed in 0 cells
summary(colMeans(counts(Mair_ds) == 0))
summary(rowMeans(counts(Mair_ds) == 0))
allzero <- rowMeans(counts(Mair_ds) == 0) == 1
Mair_ds <- Mair_ds[which(!allzero), ]

# Delete multiplets
Mair_ds <- Mair_ds[, which(colData(Mair_ds)$Sample_Tag != "Multiplet")]

# Check if we have mito genes

#getDistributionPlots()

rowData(Mair_ds)$location <- mapIds(EnsDb.Hsapiens.v86,
                                     keys=rowData(Mair_ds)$Symbol, 
                                     column="SEQNAME", keytype="GENENAME")
is.mito <- which(rowData(Mair_ds)$location=="MT")

Mair_ds = addPerCellQC(Mair_ds)

stats <- perCellQCMetrics(Mair_ds, subsets=list(Mito=which(rowData(Mair_ds)$location=="MT")))
high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
table(high.mito)

# No batch correction needed because We found that cells from different donors comingled and separated by cell type rather than by donor, suggesting that there was no cluster driven by a single donor (Figure 3A, left panel). Also, the two experimental batches (cartridges) did not show significant separation on a UMAP plot calculated using transcripts or proteins (Figure 3A, right panel).

# Discarding cells using perCellQCMetrics
df <- perCellQCMetrics(Mair_ds)
reasons <- perCellQCFilters(df)
Mair_ds <- Mair_ds[,!reasons$discard]

Mair_ds <- logNormCounts(Mair_ds)
Mair_ds <- runPCA(Mair_ds)

# Load reference data for SingleR
ref <- MonacoImmuneData()

# Perform SingleR prediction on 'sce' using the reference data
pred <- SingleR(test = Mair_ds, ref = ref, labels = ref$label.fine, assay.type.test = "logcounts")

# Assign the predicted cell types to 'sce$CellTypes'
Mair_ds$CellTypes <- factor(pred$pruned.labels)
