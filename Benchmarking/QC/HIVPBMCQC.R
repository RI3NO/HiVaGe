# Read matrix from the URL
mat <- read.table("2190_PBMC_dge.txt", header = TRUE, row.names=1)
allzero <- rowMeans(mat == 0) == 1
mat <- mat[which(!allzero), ]


# Making a plot that shows that there's clear distinction in cells that have big UMIs and low UMIs.
# Inflection is 100 UMIs so we could use default threshold.Big difference in genes expression we consider normal,
# as we have PBMC cells, most of which are specific to some function and don't express big amount of genes.
bcrank <- barcodeRanks(mat)
# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)
legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

# Deleting empty droplets using false discovery rate threshhold
e.out <- emptyDrops(mat)
summary(e.out$FDR <= 0.001)
table(colSums(mat)>100, e.out$FDR<=0.001, useNA = "ifany")
# By default all droplets with fewer than 100 UMIs are considered empty.
mat <- mat[, which(e.out$FDR <= 0.001)]

# Checking the result of filtering
bcrank <- barcodeRanks(mat)
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)
legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

# Loading the annotation dataset
library(EnsDb.Hsapiens.v86)
# Get a vector containing gene loci information
location <- mapIds(EnsDb.Hsapiens.v86,
                                    keys=rownames(mat), 
                                    column="SEQNAME", keytype="GENENAME")
# QC metric based on mito gene count
stats <- perCellQCMetrics(mat, subsets=list(Mito=which(location=="MT")))
# Find outliers
high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
# Remove them from the dataset
mat <- mat[,high.mito == FALSE]

# Get QC the other metrics
df <- perCellQCMetrics(mat)
reasons <- perCellQCFilters(df)
# Remove outliars
mat <- mat[,!reasons$discard]

# Load reference data for SingleR
ref <- MonacoImmuneData()

# Create SingleCellExperiment object based on the matrix
sce <- SingleCellExperiment(mat)
sce <- SingleCellExperiment(list(counts=mat))
sce <- logNormCounts(sce)

# Perform SingleR prediction on 'sce' using the reference data
pred <- SingleR(test = sce, ref = ref, labels = ref$label.fine, assay.type.test = "logcounts")

# Assign the predicted cell types to 'sce$CellTypes'
sce$CellTypes <- factor(pred$pruned.labels)

# Run PCA
sce <- runPCA(sce)
