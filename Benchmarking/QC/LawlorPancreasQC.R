library(scRNAseq)
library(ensembldb)
library(scater)

#loading dataset
sce.lawlor =  LawlorPancreasData()

# Making a plot that shows that there's clear distinction in cells that have big UMIs and low UMIs. Inflection is 100 UMIs so we could use default threshold. Big difference in genes expression we consider normal, as we have PBMC cells, most of which are specific to some function and don't express big amount of genes.
bcrank <- barcodeRanks(counts(sce.lawlor))

# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2)

abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)


# Delete cells that have 0 genes expressed and genes that are expressed in 0 cells
summary(colMeans(counts(sce.lawlor) == 0))
summary(rowMeans(counts(sce.lawlor) == 0))
allzero <- rowMeans(counts(sce.lawlor) == 0) == 1
sce.lawlor <- sce.lawlor[which(!allzero), ]

#annotation
library(AnnotationHub)
edb <- AnnotationHub()[["AH73881"]]
anno <- select(edb, keys=rownames(sce.lawlor), keytype="GENEID", 
               columns=c("SYMBOL", "SEQNAME"))
rowData(sce.lawlor) <- anno[match(rownames(sce.lawlor), anno[,1]),-1]

#unfiltered data
unfiltered <- sce.lawlor

#qc and filtering
stats <- perCellQCMetrics(sce.lawlor, 
                          subsets=list(Mito=which(rowData(sce.lawlor)$SEQNAME=="MT")))
qc <- quickPerCellQC(stats, percent_subsets="subsets_Mito_percent",
                     batch=sce.lawlor$`islet unos id`)
sce.lawlor <- sce.lawlor[,!qc$discard]
