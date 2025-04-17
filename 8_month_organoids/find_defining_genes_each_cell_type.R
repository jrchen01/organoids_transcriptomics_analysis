# This script aims to check the top genes defining each cluster
# April 6th, 2025

#devtools::install_github("immunogenomics/presto")
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(devtools)
library(presto)

astro_mo8 <- readRDS("/raidixshare_logg01/jchen/project/01_sarah_proj/03_seurat_files/seurat_mo8/seurat_8mo_annots_new_0306.rds")
DimPlot(astro_mo8, reduction = "umap", label = T, group.by = "cell_type") #just check
head(astro_mo8@meta.data, 3) #just check

unique(astro_mo8@meta.data$cell_type)
Idents(astro_mo8) <- "cell_type"
markers_NKCC1_pos_cells <- FindMarkers(astro_mo8, ident.1 = "NKCC1_pos_cells")
markers_ASC <- FindMarkers(astro_mo8, ident.1 = "ASC")
markers_KCC2_EX <- FindMarkers(astro_mo8, ident.1 = "KCC2_EX")
markers_GABA_INH <- FindMarkers(astro_mo8, ident.1 = "GABA_INH")
markers_OPC <- FindMarkers(astro_mo8, ident.1 = "OPC")

# use a function to extract top i gene
get_top_up_down_genes <- function(marker_df, top_n = i) {
  top_up <- marker_df %>%
    filter(avg_log2FC > 0) %>%
    arrange(p_val_adj, desc(avg_log2FC)) %>%
    slice_head(n = top_n)
  
  top_down <- marker_df %>%
    filter(avg_log2FC < 0) %>%
    arrange(p_val_adj, abs(avg_log2FC)) %>%
    slice_head(n = top_n)
  
  list(
    top_up_df = top_up,
    top_down_df = top_down,
    up_genes = rownames(top_up),
    down_genes = rownames(top_down))
}
ASC_result <- get_top_up_down_genes(markers_ASC, top_n = 10)
ASC_result$up_genes      
ASC_result$down_genes

OPC_result <- get_top_up_down_genes(markers_OPC, top_n = 10)
OPC_result$up_genes      
OPC_result$down_genes

NKCC_result <- get_top_up_down_genes(markers_NKCC1_pos_cells, top_n = 10)
NKCC_result$up_genes      
NKCC_result$down_genes

KCC2_EX_result <- get_top_up_down_genes(markers_KCC2_EX, top_n = 10)
KCC2_EX_result$up_genes      
KCC2_EX_result$down_genes

GABA_INH_result <- get_top_up_down_genes(markers_GABA_INH, top_n = 10)
GABA_INH_result$up_genes      
GABA_INH_result$down_genes
