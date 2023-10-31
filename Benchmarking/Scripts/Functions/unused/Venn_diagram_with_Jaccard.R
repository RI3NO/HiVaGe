# Load the required library for Venn diagrams
library(VennDiagram)

# Function to calculate the Jaccard index
jaccard_index <- function(set1, set2) {
  if (length(set1) == 0 || length(set2) == 0) {
    return(0)  # Return 0 if one of the sets has zero size
  }
  
  intersection <- length(intersect(set1, set2))  # Calculate the number of common elements
  union <- length(union(set1, set2))  # Calculate the number of unique elements in both sets
  
  jaccard_index <- intersection / union  # Calculate the Jaccard index
  return(jaccard_index)
}

# Read and process the BASiCS_matrix.csv file
data <- read.csv("BASiCS_matrix.csv", header = FALSE)
BASiCS_matrix <- as.matrix(data)
BASiCS_matrix <- BASiCS_matrix[2:31, ]

# Read and process the M3Drop_matrix.csv file
data <- read.csv("M3Drop_matrix.csv", header = FALSE)
M3Drop_matrix <- as.matrix(data)
M3Drop_matrix <- M3Drop_matrix[2:31, ]

# Read and process the ROGUE_matrix.csv file
data <- read.csv("ROGUE_matrix.csv", header = FALSE)
ROGUE_matrix <- as.matrix(data)
ROGUE_matrix <- ROGUE_matrix[2:31, ]

# Read and process the scran_matrix.csv file
data <- read.csv("scran_matrix.csv", header = FALSE)
scran_matrix <- as.matrix(data)
scran_matrix <- scran_matrix[2:31, ]

# Read and process the seurat_matrix.csv file
data <- read.csv("seurat_matrix.csv", header = FALSE)
seurat_matrix <- as.matrix(data)
seurat_matrix <- seurat_matrix[2:31, ]

# Define the indices of the biggest clusters to create individual Venn diagrams
biggest_clusters <- c(1, 5, 12, 13)
color_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")

# File for Jaccard index matrices
file_conn <- "Jaccard index output.txt"
if (file.exists(file_conn)) {
  file.remove(file_conn)
}
sink(file_conn, append = TRUE)

for (i in biggest_clusters) {
  methods_names = c("BASiCS", "M3Drop", "ROGUE", "scran", "seurat")
  
  # Extract gene sets for the current cluster from each matrix
  BASiCS_genes <- BASiCS_matrix[, i]
  M3Drop_genes <- M3Drop_matrix[, i]
  ROGUE_genes <- ROGUE_matrix[, i]
  scran_genes <- scran_matrix[, i]
  seurat_genes <- seurat_matrix[, i]
  
  # Create a unique filename for the Venn diagram
  constant <- "_venn_diagram.png"
  prefix <- "c"
  filename <- paste(prefix, i, constant, sep = "")
  
  # Generate a Venn diagram for the gene sets using 'venn.diagram' function
  venn.diagram(
    x = list(BASiCS_genes, M3Drop_genes, ROGUE_genes, scran_genes, seurat_genes),
    category.names = methods_names,
    filename = filename,
    output = TRUE,
    col = color_palette,   # Apply custom color palette
    fill = color_palette,  # Fill circles with custom colors
    alpha = 0.5,           # Adjust transparency of circles
    cex = 1.5              # Increase font size of labels
  )
  
  # Define the gene sets for Jaccard index matrix
  gene_sets <- list(BASiCS_genes, M3Drop_genes, ROGUE_genes, scran_genes, seurat_genes)
  num_methods <- length(gene_sets)
  
  # Create matrix with size of num_methods
  jaccard_matrix <- matrix(NA, nrow = num_methods, ncol = num_methods)
  rownames(jaccard_matrix) <- methods_names
  colnames(jaccard_matrix) <- methods_names
  
  # Calculate the Jaccard index for each pair of gene sets and fill the matrix
  for (m in 1:num_methods) {
    for (n in 1:num_methods) {
      jaccard_matrix[m, n] <- jaccard_index(gene_sets[[m]], gene_sets[[n]])
    }
  }
  # Write the Jaccard index matrix to output file
  print(paste(prefix, i, " claster:", sep=""))
  print(jaccard_matrix)
}
sink()

# Set up the color palette for the Venn diagrams
color_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")

# Create a Venn diagram for the overlapping gene sets in BASiCS, M3Drop, ROGUE, and scran
venn.diagram(
  x = list(BASiCS_matrix, M3Drop_matrix, ROGUE_matrix, scran_matrix),
  category.names = c("BASiCS", "M3Drop", "ROGUE", "scran"),
  filename = "Venn_diagram.png",
  output = TRUE,
  col = color_palette,   # Apply custom color palette
  fill = color_palette,  # Fill circles with custom colors
  alpha = 0.5,           # Adjust transparency of circles
  cex = 1.5              # Increase font size of labels
)

# Set up the color palette for the final combined Venn diagram
color_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")

# Create a combined Venn diagram for the overlapping gene sets in BASiCS, M3Drop, ROGUE, scran, and seurat
venn.diagram(
  x = list(BASiCS_matrix, M3Drop_matrix, ROGUE_matrix, scran_matrix, seurat_matrix),
  category.names = c("BASiCS", "M3Drop", "ROGUE", "scran", "seurat"),
  filename = "Whole_Venn_diagram.png",
  output = TRUE,
  col = color_palette,   # Apply custom color palette
  fill = color_palette,  # Fill circles with custom colors
  alpha = 0.5,           # Adjust transparency of circles
  cex = 1.5,             # Increase font size of labels
)
