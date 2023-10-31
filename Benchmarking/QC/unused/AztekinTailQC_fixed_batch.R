library(scRNAseq)
library(ensembldb)
library(scater)
library(xlaevis.db)
library(org.Xl.eg.db)
library(AnnotationHub)
library(scuttle)
library(DropletUtils)
library(scDblFinder)

#loading dataset
Aztekin_ds <- AztekinTailData()
# No need for batch correction because there's clustering 
bcrank <- barcodeRanks(counts(Aztekin_ds))
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2,
     ylim=c(1, max(bcrank$total[uniq])))

abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

# Aztekin_ds %>% 
#   counts() %>%
#   class()
# # Deleting empty droplets using false discovery rate threshhold
# e.out <- emptyDrops(counts(Aztekin_ds))
# summary(e.out$FDR <= 0.001)
# table(colSums(counts(Aztekin_ds))>100, e.out$FDR<=0.001, useNA = "ifany")
# # By default all droplets with fewer than 100 UMIs are considered empty. But there's no droplets fewer than 1000 UMI.
# Aztekin_ds <- Aztekin_ds[, which(e.out$FDR <= 0.001)]

# Delete cells that have 0 genes expressed and genes that are expressed in 0 cells
# summary(colMeans(counts(Aztekin_ds) == 0))
# summary(rowMeans(counts(Aztekin_ds) == 0))
allzero <- rowMeans(counts(Aztekin_ds) == 0) == 1
Aztekin_ds <- Aztekin_ds[which(!allzero), ]


Aztekin_ds <- logNormCounts(Aztekin_ds)

# Dimensionality reduction with PCA
# Run PCA on 'sce'
Aztekin_ds <- runPCA(Aztekin_ds)
# Compute doublet density using the reduced dimensions from PCA
dbl.dens <- computeDoubletDensity(Aztekin_ds, d = ncol(reducedDim(Aztekin_ds)))
# Add doublet scores to the 'sce' object
Aztekin_ds$DoubletScore <- dbl.dens
# Remove doublets from 'sce' based on the identified doublet calls
# Identify doublets using doubletThresholding function
dbl.calls <- doubletThresholding(data.frame(score = dbl.dens),
                                 method = "griffiths", returnType = "call")
# table(dbl.calls == "singlet")
# Remove doublets from 'sce' based on the identified doublet calls
Aztekin_ds <- Aztekin_ds[, dbl.calls == "singlet"]
# Perform PCA again
Aztekin_ds <- runPCA(Aztekin_ds)



ah <- AnnotationHub()
annotation = ah[['AH111581']]
annotated_dataset <- select(annotation, keys = rownames(Aztekin_ds), keytype = "SYMBOL", column = 'GENENAME')
is.mito <- grepl("mitoch", annotated_dataset$GENENAME, ignore.case = TRUE)

Aztekin_ds <- addPerCellQC(Aztekin_ds)
Aztekin_ds <- addPerCellQC(Aztekin_ds, subsets=list(Mito=is.mito))

batch_ = paste0(Aztekin_ds$Condition, "_bat-", Aztekin_ds$batch)
Aztekin_ds$new_batch = batch_

discard.ercc <- isOutlier(Aztekin_ds$subsets_Mito_percent,
                          type="higher", batch=batch_)
ercc.thresholds <- attr(discard.ercc, "thresholds")["higher",]
ercc.thresholds
names(ercc.thresholds)[isOutlier(ercc.thresholds, type="higher")]
discard.ercc2 <- isOutlier(Aztekin_ds$subsets_Mito_percent,
                           type="higher", batch=batch_,
                           subset=!(batch_ %in% names(ercc.thresholds)[isOutlier(ercc.thresholds, type="higher")]))
batch.reasons <- perCellQCFilters(Aztekin_ds, batch=batch_,
                                  sub.fields=c("subsets_Mito_percent"))

Aztekin_ds = Aztekin_ds[,!batch.reasons$discard]


df <- perCellQCMetrics(Aztekin_ds)
reasons <- perCellQCFilters(df)
Aztekin_ds <- Aztekin_ds[,!reasons$discard]


Aztekin_ds <- logNormCounts(Aztekin_ds)
Aztekin_ds <- runPCA(Aztekin_ds)
