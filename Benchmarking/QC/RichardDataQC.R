library(scRNAseq)
Richard_ds <- RichardTCellData(location = TRUE)

# Making a plot that shows that there's clear distinction in cells that have big UMIs and low UMIs. Inflection is 100 UMIs so we could use default threshold. Big difference in genes expression we consider normal, as we have PBMC cells, most of which are specific to some function and don't express big amount of genes.
bcrank <- barcodeRanks(counts(Richard_ds))

# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2)

abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)


allzero <- rowMeans(counts(Richard_ds) == 0) == 1
Richard_ds <- Richard_ds[which(!allzero), ]

#cleaning up in accordance with metadata
Richard_ds = Richard_ds[,Richard_ds$`single cell quality` == "OK"]

Richard_ds <- logNormCounts(Richard_ds)
Richard_ds <- runPCA(Richard_ds)

# Load reference data for SingleR
ref <- ImmGenData(ensembl = TRUE)

# Perform SingleR prediction on 'Richard_ds' using the reference data
pred <- SingleR(test = Richard_ds, ref = ref, labels = ref$label.fine, assay.type.test = "logcounts")

# Assign the predicted cell types to 'Richard_ds$CellTypes'
Richard_ds$CellTypes <- factor(pred$pruned.labels)

Richard_ds$CellTypes

Richard_ds$new_batch = paste0(Richard_ds$time, Richard_ds$stimulus)