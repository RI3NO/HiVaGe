All files that were generated during evaluation of 5 HVGs identification tools (seurat, scran, BASiCS, M3Drop, ROGUE) of Mair dataset from scRNAseq R package

clusters - Louvain clusters that were generated before HVG identification

method_matrix.csv - matrices where columns are different Louvain clusters and contain HVGs (where method is a HVG identification method that was used to generate HVGs in each cluster)

c*_venn_diagram - Venn diagrams from HVGs that were obtained from 4 biggest Louvain clusters (where * is a number of cluster) 

cell_percentages.csv - table in which 1-st row contains cell type names, 2-nd contains cell type precentages, obtained using all genes and 3-rd contains cell type precentages that were obtained, using most convenient (derived from Venn diagrams) method, seurat. Cell types were obtained by using MonacoImmuneData() from SingleR R package as external reference

Presentation.pdf - presentation, that was generated during EEBG school, showing this work 
 
tSNE_clusters.png - tSNE plot of Louvain clusters that were generated before HVG identification

Venn_diagram.R - code that imports method_matrix.csv and outports Whole_Venn_diagram.png

Whole_Venn_diagram.png - Venn diagram from HVGs that were obtained all Louvain clusters 
