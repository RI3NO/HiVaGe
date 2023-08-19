library(scRNAseq)
library(ensembldb)
library(scater)
library(scuttle)
library(DropletUtils)
library(AnnotationHub)

# Load MairPBMCData dataset and store it in 'sce' variable
sce<- MairPBMCData(ensembl = TRUE)

sce

#annotation
ah <- AnnotationHub()

annotation <- ah[["AH109606"]]


location <- mapIds(annotation, keys= rownames(sce), column = 'SEQNAME', keytype="GENEID")



# Filtering
# Remove empty droplets from the dataset
set.seed(42)
e.out <- emptyDrops(counts(sce))
summary(e.out$FDR <= 0.001)

sce <- sce[,which(e.out$FDR <= 0.001)]


unfiltered <- sce




#quality control 

stats <- perCellQCMetrics(sce, subsets=list(Mito=which(location=="MT")))
stats


high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
low.adt <- stats$`altexps_adt_detected` < nrow(altExp(sce))/2

discard <- high.mito | low.adt
sce <- sce[,!discard]





source("HiVaGe.R")

# Call the function
test <- HiVaGe(sce_filtered, flavour = "BASiCS", percentile = 0.5)
print(test)




# #removing cells with 0 counts
# 
# # Check which genes have zero counts in all cells
# genes_with_zero_counts <- rowSums(counts(sce)) == 0
# 
# # Create a new expression matrix excluding the zero-count genes
# sce <- sce[!genes_with_zero_counts, ]


#Remove transcripts with low total counts
total_counts_threshold <- 10  # Set your desired total counts threshold here, originally was 10

#Calculate total counts for each transcript
transcript_total_counts <- rowSums(counts(sce))

# Subset transcripts with total counts above the threshold
sce <- sce[transcript_total_counts >= total_counts_threshold, ]

# # Step 2: Remove transcripts expressed in a few cells
# expressed_in_few_cells_threshold <- 2  # Set your desired expressed in few cells threshold here (was 5)
# 
# # Calculate the number of cells where each transcript is expressed
# transcript_expression_cells <- rowSums(counts(sce_filtered_counts) > 0)
# 
# # Subset transcripts expressed in more cells than the threshold
# sce_filtered <- sce_filtered_counts[transcript_expression_cells >= expressed_in_few_cells_threshold, ]

# Optional: Update other assay data (e.g., log-normalized counts or normalized counts) in the filtered sce object if needed
# Example: sce_filtered$logcounts <- log(counts(sce_filtered) + 1)


#creating graphs


# summary(low.adt)
# summary(discard)
#           
# colData(unfiltered) <- cbind(colData(unfiltered), stats)
# unfiltered$discard <- discard
# 
# gridExtra::grid.arrange(
#   plotColData(unfiltered, y="sum", colour_by="discard") +
#     scale_y_log10() + ggtitle("Total count"),
#   plotColData(unfiltered, y="detected", colour_by="discard") +
#     scale_y_log10() + ggtitle("Detected features"),
#   plotColData(unfiltered, y="subsets_Mito_percent",
#               colour_by="discard") + ggtitle("Mito percent"),
#   plotColData(unfiltered, y="altexps_adt_detected",
#               colour_by="discard") + ggtitle("ADT detected"),
#   ncol=2
# )


# #export dataset 
# 
# # Access the gene expression data from the SCE object
# expression_data <- assays(sce)$counts
# 
# # Convert the gene expression data into a data frame
# expression_df <- as.data.frame(expression_data)
# 
# # Export the data frame as a CSV file
# write.csv(expression_df, file = "C:/Users/Asyandel/Documents/r_projects/HVG/mairpbmc_qcdone.csv", row.names = FALSE)
# 

