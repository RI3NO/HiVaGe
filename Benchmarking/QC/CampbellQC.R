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
Campbell_ds = CampbellBrainData()

rowData(Campbell_ds)$is_spike = grepl("ERCC", rownames(Campbell_ds))
altExp(Campbell_ds, "ERCC") = Campbell_ds[rowData(Campbell_ds)$is_spike]
Campbell_ds = Campbell_ds[!rowData(Campbell_ds)$is_spike]
# No need for batch correction because there's clustering 
bcrank <- barcodeRanks(counts(Campbell_ds))
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2,
     ylim=c(1, max(bcrank$total[uniq])))

abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

# Campbell_ds %>% 
#   counts() %>%
#   class()
# # Deleting empty droplets using false discovery rate threshhold
# e.out <- emptyDrops(counts(Campbell_ds))
# summary(e.out$FDR <= 0.001)
# table(colSums(counts(Campbell_ds))>100, e.out$FDR<=0.001, useNA = "ifany")
# # By default all droplets with fewer than 100 UMIs are considered empty. But there's no droplets fewer than 1000 UMI.
# Campbell_ds <- Campbell_ds[, which(e.out$FDR <= 0.001)]

allzero <- rowMeans(counts(Campbell_ds) == 0) == 1
# Create a table showing the number of cells with all-zero counts
table(allzero)
# Remove cells with all-zero counts for all genes from 'sce'
Campbell_ds <- Campbell_ds[which(!allzero),]


Campbell_ds <- logNormCounts(Campbell_ds)

# Dimensionality reduction with PCA
# Run PCA on 'sce'
Campbell_ds <- runPCA(Campbell_ds)
# Compute doublet density using the reduced dimensions from PCA
dbl.dens <- computeDoubletDensity(Campbell_ds, d = ncol(reducedDim(Campbell_ds)))
# Add doublet scores to the 'sce' object
Campbell_ds$DoubletScore <- dbl.dens
# Remove doublets from 'sce' based on the identified doublet calls
# Identify doublets using doubletThresholding function
dbl.calls <- doubletThresholding(data.frame(score = dbl.dens),
                                 method = "griffiths", returnType = "call")
# table(dbl.calls == "singlet")
# Remove doublets from 'sce' based on the identified doublet calls
Campbell_ds <- Campbell_ds[, dbl.calls == "singlet"]
# Perform PCA again
Campbell_ds <- runPCA(Campbell_ds)

#annotation
mm_87 = AnnotationHub()[["AH53222"]]
chr.loc = mapIds(mm_87,
       keys=rownames(rowData(Campbell_ds)), 
       column="SEQNAME", keytype="SYMBOL")
is.mito <- which(chr.loc=="MT")

Campbell_ds <- addPerCellQC(Campbell_ds)
Campbell_ds <- addPerCellQC(Campbell_ds, subsets=list(Mito=is.mito))

discard.ercc <- isOutlier(Campbell_ds$subsets_Mito_percent,
                          type="higher", batch=Campbell_ds$group)
ercc.thresholds <- attr(discard.ercc, "thresholds")["higher",]
ercc.thresholds
names(ercc.thresholds)[isOutlier(ercc.thresholds, type="higher")]
discard.ercc2 <- isOutlier(Campbell_ds$subsets_Mito_percent,
                           type="higher", batch=Campbell_ds$group,
                           subset=!(Campbell_ds$group %in% names(ercc.thresholds)[isOutlier(ercc.thresholds, type="higher")]))
batch.reasons <- perCellQCFilters(Campbell_ds, batch=Campbell_ds$group,
                                  sub.fields=c("subsets_Mito_percent"))

Campbell_ds = Campbell_ds[,!batch.reasons$discard]


df <- perCellQCMetrics(Campbell_ds)
reasons <- perCellQCFilters(df)
Campbell_ds <- Campbell_ds[,!reasons$discard]


Campbell_ds <- logNormCounts(Campbell_ds)
Campbell_ds <- runPCA(Campbell_ds)

new_rownames = gsub("_", "-", rownames(Campbell_ds))
rownames(Campbell_ds) = make.names(new_rownames, unique=TRUE)
