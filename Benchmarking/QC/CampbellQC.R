library(scRNAseq)
library(ensembldb)
library(scater)
library(xlaevis.db)
library(org.Xl.eg.db)
library(AnnotationHub)
library(scuttle)
library(DropletUtils)
library(scDblFinder)

# Loading the dataset from scRNAseq library
Campbell_ds = CampbellBrainData()

# Extract ERCC data from counts table
rowData(Campbell_ds)$is_spike = grepl("ERCC", rownames(Campbell_ds))
# And put it into altExp field of SingleCellExperiment
altExp(Campbell_ds, "ERCC") = Campbell_ds[rowData(Campbell_ds)$is_spike]
# Remove ERCC data from counts table
Campbell_ds = Campbell_ds[!rowData(Campbell_ds)$is_spike]

# No need for batch correction because it was already done by the authors

## Knee plot preparation
bcrank <- barcodeRanks(counts(Campbell_ds))
# Only showing unique points for plotting speed
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2,
     ylim=c(1, max(bcrank$total[uniq])))
abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)
legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

# Removing empty drops is not required as per how the knee plot looks
# In addition, e.out function itself tells this dataset has already been filtered for empty droplets

allzero <- rowMeans(counts(Campbell_ds) == 0) == 1
# Create a table showing the number of cells with all-zero counts
table(allzero)
# Remove cells with all-zero counts for all genes from 'sce'
Campbell_ds <- Campbell_ds[which(!allzero),]

# Logcounts for PCA
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

# Loading the annotation dataset
mm_87 = AnnotationHub()[["AH53222"]]
# Get a vector containing gene loci information
chr.loc = mapIds(mm_87,
       keys=rownames(rowData(Campbell_ds)), 
       column="SEQNAME", keytype="SYMBOL")
# Vector containing information about if gene is mitochondrial
is.mito <- which(chr.loc=="MT")

# Get QC metrics for Campbell_ds
Campbell_ds <- addPerCellQC(Campbell_ds)
# Get QC metric based on mito genes
Campbell_ds <- addPerCellQC(Campbell_ds, subsets=list(Mito=is.mito))

# Perform outlier detection based on mito gene count in different batches
discard.ercc <- isOutlier(Campbell_ds$subsets_Mito_percent,
                          type="higher", batch=Campbell_ds$group)
# Get the thresholds into a separate value
ercc.thresholds <- attr(discard.ercc, "thresholds")["higher",]
# Get names of outlier batches
names(ercc.thresholds)[isOutlier(ercc.thresholds, type="higher")]
# No outlier batches

# Get outliers based on all batches 
batch.reasons <- perCellQCFilters(Campbell_ds, batch=Campbell_ds$group,
                                  sub.fields=c("subsets_Mito_percent"))
# Remove outliers from the dataset
Campbell_ds = Campbell_ds[,!batch.reasons$discard]

# Get QC the other metrics
df <- perCellQCMetrics(Campbell_ds)
reasons <- perCellQCFilters(df)
# Remove outliars
Campbell_ds <- Campbell_ds[,!reasons$discard]

# Fix rownames (gene names)
new_rownames = gsub("_", "-", rownames(Campbell_ds))
rownames(Campbell_ds) = make.names(new_rownames, unique=TRUE)
