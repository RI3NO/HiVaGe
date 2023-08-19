library(scRNAseq)
library(ensembldb)
library(scater)
library(scuttle)
library(DropletUtils)
library(EnsDb.Hsapiens.v86)

#loading dataset
sce.stoe <- StoeckiusHashingData(type = 'pbmc', mode= 'human', ensembl = TRUE)

sce.stoe



#removing empty droplets(or with high discovery rate)
set.seed(100)
e.out <- emptyDrops(counts(sce.stoe))
sce.stoe <- sce.stoe[,which(e.out$FDR <= 0.001)]

#creating a copy
unfiltered <- sce.stoe

#maping genes to chromosomes?)

ah <- AnnotationHub()

annotation <- ah[["AH109606"]]

location_geneids <- mapIds(annotation, keys= rownames(sce.stoe), column = 'SEQNAME', keytype="GENEID")

# location_genenames <- mapIds(EnsDb.Hsapiens.v86, keys= rowData(sce.stoe)$originalName,
#                               column="SEQNAME", keytype="GENENAME")

#actually doing QC on mito genes
stats <- perCellQCMetrics(sce.stoe, subsets=list(Mito=which(location_geneids == "MT")))
high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")

#removing NAs
high.mito <- high.mito[complete.cases(high.mito)]

sce.stoe <- sce.stoe[,!high.mito]


#Remove transcripts with low total counts
total_counts_threshold <- 10  # Set your desired total counts threshold here, originally was 10

#Calculate total counts for each transcript
transcript_total_counts <- rowSums(counts(sce.stoe))

# Subset transcripts with total counts above the threshold
sce <- sce.stoe[transcript_total_counts >= total_counts_threshold, ]

# summary(high.mito)
# 
# 
# 
# colData(unfiltered) <- cbind(colData(unfiltered), stats)
# unfiltered$discard <- high.mito
# 
# gridExtra::grid.arrange(
#   plotColData(unfiltered, y="sum", colour_by="discard") +
#     scale_y_log10() + ggtitle("Total count"),
#   plotColData(unfiltered, y="detected", colour_by="discard") +
#     scale_y_log10() + ggtitle("Detected features"),
#   plotColData(unfiltered, y="subsets_Mito_percent",
#               colour_by="discard") + ggtitle("Mito percent"),
#   ncol=2
#)
