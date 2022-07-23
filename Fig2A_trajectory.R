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

seurat_data <- Read10X("/bigdata/godziklab/shared/Xinru/302011/Fig2A_Outs_Folder/filtered_feature_bc_matrix/")
fig2a <- CreateSeuratObject(counts = seurat_data)
barcode_cluster <- read.csv("/bigdata/godziklab/shared/Xinru/302011/Fig2A_Outs_Folder/LibraryID.csv")

### add cluster info to the single cell data
meta <- fig2a@meta.data
meta$Barcode <- row.names(meta)
meta2 <- left_join(meta, barcode_cluster, by = "Barcode")

dim(meta)
dim(meta2)
row.names(meta2) <- meta2$Barcode

meta2$Genotype <- ifelse(grepl("mGluR\\-", meta2$LibraryID), "mGluR-", "")
meta2$Genotype <- ifelse(grepl("mGluR\\+", meta2$LibraryID), "mGluR+", meta2$Genotype)
table(meta2$Genotype)

# sample1 sample2 sample3
# 283     279     552
# 1114 cells
fig2a@meta.data <- meta2
table(fig2a@meta.data$Genotype)

# Trajectory analysis -----------------------------------------------------

exp_mtx <- as.matrix(GetAssayData(fig2a[["RNA"]], slot = "counts"))
cell_meta <- fig2a@meta.data
gene_ano <- data.frame(gene_short_name=row.names(exp_mtx))
rownames(gene_ano) <- gene_ano$gene_short_name
cds <- new_cell_data_set(exp_mtx,
                         cell_metadata = cell_meta,
                         gene_metadata = gene_ano)

cds <- preprocess_cds(cds, num_dim = 100)
cds <- align_cds(cds, alignment_group = "LibraryID")
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

cds@colData$Cluster <- factor(cds@colData$Genotype,
                              levels= c('mGluR-', 'mGluR+'))
setwd("/bigdata/godziklab/shared/Xinru/302011/Monocle3_output")
dpi = 300
png(file = "Monocle3_Traject_2A_v1.png", width = dpi * 9,height = dpi * 6,units = "px",res = dpi,type = 'cairo')
plot_cells(cds,
           color_cells_by = "Genotype",
           cell_size = 0.5,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=1.5) +
  scale_color_manual(values = c('#4A90E2', '#F5A623'))
dev.off()


setwd("/bigdata/godziklab/shared/Xinru/302011/Monocle3_output")
dpi = 300
png(file = "Monocle3_Traject_2A_v2.png", width = dpi * 6,height = dpi * 4,units = "px",res = dpi,type = 'cairo')
plot_cells(cds,
           color_cells_by = "Genotype",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           label_principal_points = TRUE,
           graph_label_size=1.5) 
dev.off()


cds <- order_cells(cds, 
                   reduction_method = "UMAP",
                   root_pr_nodes = "Y_1")

dpi = 300
png(file = "Monocle3_Traject_2A_v3.png", width = dpi * 6,height = dpi * 4,units = "px",res = dpi,type = 'cairo')
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5) 
dev.off()
