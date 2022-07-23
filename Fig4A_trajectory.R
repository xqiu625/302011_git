.libPaths( c( "/bigdata/godziklab/shared/Xinru/R" , .libPaths() ) )
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(stringr)
library(dplyr)
library(plyr)
library(monocle3)

seurat_data <- Read10X("/bigdata/godziklab/shared/Xinru/302011/Fig4A_Outs_Folder/filtered_feature_bc_matrix/")
fig4a <- CreateSeuratObject(counts = seurat_data)
barcode_cluster <- read.csv("/bigdata/godziklab/shared/Xinru/302011/Fig4A_Outs_Folder/K-Means\ 4.csv")


### add cluster info to the single cell data
meta <- fig4a@meta.data
meta$Barcode <- row.names(meta)
meta2 <- left_join(meta, barcode_cluster, by = "Barcode")

dim(meta2)
row.names(meta2) <- meta2$Barcode
colnames(meta2)[colnames(meta2)=="X4"] <- "Cluster"
meta2$Sample.ID <- ifelse(grepl("3$", meta2$Barcode), "sample3", "")
meta2$Sample.ID <- ifelse(grepl("2$", meta2$Barcode), "sample2", meta2$Sample.ID)
meta2$Sample.ID <- ifelse(grepl("1$", meta2$Barcode), "sample1", meta2$Sample.ID)
table(meta2$Sample.ID)

# sample1 sample2 sample3
# 324     334     570
fig4a@meta.data <- meta2
table(fig4a@meta.data$Cluster)

fig4a_filtered <- subset(fig4a, Cluster %in% c('Cluster 1', 'Cluster 2', 'Cluster 3'))
# Trajectory analysis -----------------------------------------------------

exp_mtx <- as.matrix(GetAssayData(fig4a_filtered[["RNA"]], slot = "counts"))
cell_meta <- fig4a_filtered@meta.data
gene_ano <- data.frame(gene_short_name=row.names(exp_mtx))
rownames(gene_ano) <- gene_ano$gene_short_name
cds <- new_cell_data_set(exp_mtx,
                         cell_metadata = cell_meta,
                         gene_metadata = gene_ano)

cds <- preprocess_cds(cds, num_dim = 100)
cds <- align_cds(cds, alignment_group = "Sample.ID")
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

cds@colData$Cluster <- factor(cds@colData$Cluster,
                              levels= c('Cluster 1', 'Cluster 2', 'Cluster 3'))
setwd("/bigdata/godziklab/shared/Xinru/302011/Monocle3_output")
dpi = 300
png(file = "Monocle3_Traject_4A_v1.png", width = dpi * 9,height = dpi * 6,units = "px",res = dpi,type = 'cairo')
plot_cells(cds,
           color_cells_by = "Cluster",
           cell_size = 0.5,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=1.5) +
  scale_color_manual(values = c('#9013FE', '#7ED321', '#D0021B'))
dev.off()


setwd("/bigdata/godziklab/shared/Xinru/302011/Monocle3_output")
dpi = 300
png(file = "Monocle3_Traject_4A_v2.png", width = dpi * 6,height = dpi * 4,units = "px",res = dpi,type = 'cairo')
plot_cells(cds,
           color_cells_by = "Cluster",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           label_principal_points = TRUE,
           graph_label_size=1.5) 
dev.off()


cds <- order_cells(cds, 
                   reduction_method = "UMAP",
                   root_pr_nodes = "Y_3")

setwd("/bigdata/godziklab/shared/Xinru/302011/Monocle3_output")
dpi = 300
png(file = "Monocle3_Traject_4A_v3.png", width = dpi * 6,height = dpi * 4,units = "px",res = dpi,type = 'cairo')
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5) 
dev.off()
