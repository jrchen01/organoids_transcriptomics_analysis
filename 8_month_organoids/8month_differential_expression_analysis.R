#!/usr/bin/env Rscript

# Differential expression analysis using DESeq2
## Pay attention to the distribution of gender

library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)
library(DESeq2)
library(tidyverse)
library(ashr)

astro_mo8 <- readRDS("/raidixshare_logg01/jchen/project/01_sarah_proj/03_seurat_files/seurat_mo8/seurat_8mo_annots_0411.rds")
DimPlot(astro_mo8, reduction = "umap", label = T, group.by = "cell_type") #just check
seurat <- astro_mo8

prepare_dds <- function(seurat, cell_type) { # Prepare DESeq2 dds object function
  p <- seurat@meta.data
  p <- p[p$cell_type %in% c(cell_type), ] #select cells from only one of your cell types
  ids <- table(p$patient_id) # donor is the column name with the cell-line ID
  ids <- names(ids[ids >= 10]) # Select donors with at least 10 cells
  p <- p[p$patient_id %in% ids, ]
  counts <- seurat@assays$RNA$counts #raw counts
  counts <- counts[, rownames(p)]
  sce <- SingleCellExperiment(assays = list(counts = counts),
                              colData = p[, c("patient_id"), drop = FALSE])
  groups <- SingleCellExperiment::colData(sce)
  pb <- Matrix.utils::aggregate.Matrix(t(counts(sce)), 
                                       groupings = groups, fun = "sum") # Aggregate across cluster-sample groups
  df <- as.data.frame(t(pb))
  p <- p[!duplicated(p$patient_id), ] # Remove duplicates so that p contains only a unique row for each donor (metadata per sample, not per cell anymore)
  rownames(p) <- p$patient_id
  p$diagnosis <- factor(p$diagnosis)
  p$diagnosis <- relevel(p$diagnosis, ref = "NC") # reference for log2 Fold-change
  p <- p[colnames(df), ]
  # Checking data consistency
  if (!identical(rownames(p), colnames(df))) {
    stop("Error: rownames of p and colnames of df do not match!")
  }
  if (nrow(p) != ncol(df)) {
    stop("Error: Number of rows in p does not match number of columns in df!")
  }
  #print(table(p$diagnosis))
  return(list(df = df, p = p, diagnosis = p$diagnosis))
}

unique(seurat@meta.data$cell_type) #Cell types
EX_result <- prepare_dds(seurat, "KCC2_EX") # KCC2+ excitatory neurons
print(table(EX_result$diagnosis))
table(EX_result$p$diagnosis, EX_result$p$sex) #check the distribution of sex - imbalanced
dds_EX <- DESeqDataSetFromMatrix(EX_result$df, 
                                 colData = EX_result$p, 
                                 design = ~ diagnosis + sex)
smallestGroupSize <- 2 
keep_EX <- rowSums(counts(dds_EX) >= 10) >= smallestGroupSize 
dds_EX <- dds_EX[keep_EX,] # keep genes with at least 10 counts in two samples
dds_EX <- DESeq(dds_EX)
resultsNames(dds_EX) # just check
res_EX <- as.data.frame(results(dds_EX, name = "diagnosis_AD_vs_NC"))
# Turn the results object into a tibble for use with tidyverse functions
res_tbl_EX <- res_EX %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()
res_tbl_EX <- res_tbl_EX %>% arrange(padj)
head(res_tbl_EX, 10)
(sig_res_EX <- dplyr::filter(res_tbl_EX, padj < 0.05) %>% # any significant results
  dplyr::arrange(padj))
gene_list_EX <- res_tbl_EX %>% # save gene list (p-val < 0.05)
  filter(pvalue < 0.05) %>%
  pull(gene) 
writeLines(gene_list_EX, "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/Seurat_8mo/02_pseudobulk_de/deseq2_de_list/gene_list_EX.txt")

NKCC1_result <- prepare_dds(seurat, "NKCC1_pos_cells") # NKCC1 positive cells / NKCC1+ neurons
print(table(NKCC1_result$diagnosis))
table(NKCC1_result$p$diagnosis, NKCC1_result$p$sex) # imbalanced
dds_NKCC1 <- DESeqDataSetFromMatrix(NKCC1_result$df, colData = NKCC1_result$p, design = ~ diagnosis + sex)
smallestGroupSize <- 2 
keep_NKCC1 <- rowSums(counts(dds_NKCC1) >= 10) >= smallestGroupSize 
dds_NKCC1 <- dds_NKCC1[keep_NKCC1,] # keep genes with at least 10 counts in two samples
dds_NKCC1 <- DESeq(dds_NKCC1)
res_NKCC1 <- as.data.frame(results(dds_NKCC1, name = "diagnosis_AD_vs_NC"))
res_tbl_NKCC1 <- res_NKCC1 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()
res_tbl_NKCC1 <- res_tbl_NKCC1 %>% arrange(pvalue)
head(res_tbl_NKCC1)
(sig_res_NKCC1 <- dplyr::filter(res_tbl_NKCC1, padj < 0.05) %>% #Subset the significant results
  dplyr::arrange(padj)) # none
gene_list_NKCC1 <- res_tbl_NKCC1 %>% #Save gene list (p-val < 0.05)
  filter(pvalue < 0.05) %>%
  pull(gene) 
writeLines(gene_list_NKCC1, "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/Seurat_8mo/02_pseudobulk_de/deseq2_de_list/gene_list_NKCC1.txt")

ASC_result <- prepare_dds(seurat, "ASC") # Astrocytes
print(table(ASC_result$diagnosis))
table(ASC_result$p$diagnosis, ASC_result$p$sex) # imbalanced
dds_ASC <- DESeqDataSetFromMatrix(ASC_result$df, colData = ASC_result$p, design = ~ diagnosis + sex)
smallestGroupSize <- 2 
keep_ASC <- rowSums(counts(dds_ASC) >= 10) >= smallestGroupSize 
dds_ASC <- dds_ASC[keep_ASC,] # keep genes with at least 10 counts in two samples
dds_ASC <- DESeq(dds_ASC)
res_ASC <- as.data.frame(results(dds_ASC, name = "diagnosis_AD_vs_NC"))
res_tbl_ASC <- res_ASC %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()
res_tbl_ASC <- res_tbl_ASC %>% arrange(pvalue)
head(res_tbl_ASC)
(sig_res_ASC <- dplyr::filter(res_tbl_ASC, padj < 0.05) %>%
  dplyr::arrange(padj)) # none
gene_list_ASC <- res_tbl_ASC %>% # save gene list (p-val < 0.05)
  filter(pvalue < 0.05) %>%
  pull(gene) 
writeLines(gene_list_ASC, "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/Seurat_8mo/02_pseudobulk_de/deseq2_de_list/gene_list_ASC.txt")

OPC_result <- prepare_dds(seurat, "OPC_ODC") # Oligodendrocytes and oligodendrocyte precursor cells 
print(table(OPC_result$diagnosis))
table(OPC_result$p$diagnosis, OPC_result$p$sex) # imbalanced
dds_OPC <- DESeqDataSetFromMatrix(OPC_result$df, colData = OPC_result$p, design = ~ diagnosis + sex)
smallestGroupSize <- 2 
keep_OPC <- rowSums(counts(dds_OPC) >= 10) >= smallestGroupSize 
dds_OPC <- dds_OPC[keep_OPC,] # keep genes with at least 10 counts in two samples
dds_OPC <- DESeq(dds_OPC)
res_OPC <- as.data.frame(results(dds_OPC, name = "diagnosis_AD_vs_NC"))
res_tbl_OPC <- res_OPC %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()
res_tbl_OPC <- res_tbl_OPC %>% arrange(pvalue)
head(res_tbl_OPC)
(sig_res_OPC <- dplyr::filter(res_tbl_OPC, padj < 0.05) %>%
  dplyr::arrange(padj)) # none
gene_list_OPC <- res_tbl_OPC %>% # save gene list (p-val < 0.05)
  filter(pvalue < 0.05) %>%
  pull(gene) 
writeLines(gene_list_OPC, "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/Seurat_8mo/02_pseudobulk_de/deseq2_de_list/gene_list_OPC.txt")

INH_result <- prepare_dds(seurat, "GABA_INH") # GABAergic inhibitory neurons
print(table(INH_result$diagnosis))
table(INH_result$p$diagnosis, INH_result$p$sex) # imbalanced
dds_INH <- DESeqDataSetFromMatrix(INH_result$df, colData = INH_result$p, design = ~ diagnosis + sex)
smallestGroupSize <- 2 
keep_INH <- rowSums(counts(dds_INH) >= 10) >= smallestGroupSize 
dds_INH <- dds_INH[keep_INH,] # keep genes with at least 10 counts in two samples
dds_INH <- DESeq(dds_INH)
res_INH <- as.data.frame(results(dds_INH, name = "diagnosis_AD_vs_NC"))
res_tbl_INH <- res_INH %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()
res_tbl_INH <- res_tbl_INH %>% arrange(pvalue)
head(res_tbl_INH)
(sig_res_INH <- dplyr::filter(res_tbl_INH, padj < 0.05) %>%
    dplyr::arrange(padj)) # none
gene_list_INH <- res_tbl_INH %>% # save gene list (p-val < 0.05)
  filter(pvalue < 0.05) %>%
  pull(gene) 
writeLines(gene_list_INH, "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/Seurat_8mo/02_pseudobulk_de/deseq2_de_list/gene_list_INH.txt")
