library(scRNAseq)
library(ensembldb)
library(scater)
library(xlaevis.db)
library(org.Xl.eg.db)
library(AnnotationHub)
library(scuttle)
library(DropletUtils)

# Loading the dataset from scRNAseq library
Aztekin_ds <- AztekinTailData()
# Batch correction was already performed for this data 

# Knee plot
bcrank <- barcodeRanks(counts(Aztekin_ds))
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2,
     ylim=c(1, max(bcrank$total[uniq])))
abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)
legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

# Deleting empty droplets using false discovery rate threshhold
e.out <- emptyDrops(counts(Aztekin_ds))
summary(e.out$FDR <= 0.001)
table(colSums(counts(Aztekin_ds))>100, e.out$FDR<=0.001, useNA = "ifany")
# By default all droplets with fewer than 100 UMIs are considered empty. But there's no droplets fewer than 1000 UMI.
Aztekin_ds <- Aztekin_ds[, which(e.out$FDR <= 0.001)]

# Delete cells that have 0 genes expressed and genes that are expressed in 0 cells
summary(colMeans(counts(Aztekin_ds) == 0))
summary(rowMeans(counts(Aztekin_ds) == 0))
allzero <- rowMeans(counts(Aztekin_ds) == 0) == 1
Aztekin_ds <- Aztekin_ds[which(!allzero), ]

## Annotate genes for future QC using mitochondrial genes
ah <- AnnotationHub()
# Load the annotation for frog cells
annotation = ah[['AH111581']]
annotated_dataset <- select(annotation, keys = rownames(sce.aztekin), keytype = "SYMBOL", column = 'GENENAME')
# Vector containing information about if gene is mitochondrial
is.mito <- grepl("mitoch", annotated_dataset$GENENAME, ignore.case = TRUE)

# Get QC metrics based on mito genes
stats <- perCellQCMetrics(sce.aztekin, subsets=list(Mito=is.mito))
# Cells with high mito are outliers
high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
# Remove said cells
sce.aztekin <- sce.aztekin[,!high.mito]

# Get QC the other metrics
df <- perCellQCMetrics(sce.aztekin)
reasons <- perCellQCFilters(df)
# Remove outliars
sce.aztekin <- sce.aztekin[,!reasons$discard]




