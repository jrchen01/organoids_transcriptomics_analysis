#!/usr/bin/env Rscript

# Previously we saw a downregulation in kccn1 and maybe kccn2 (sk channels) in AD GAD1+ cell population.
# Look and see if the inhibitory neurons or maybe gad1+ AD cells still look like they have this downregulation with the new clustering.
# Use DESeq2 and pesudobulk for each patient/cell line for comparison
# April 6, 2025

library(Seurat)
library(dplyr)
library(ggplot2)

astro_mo8 <- readRDS("/raidixshare_logg01/jchen/project/01_sarah_proj/03_seurat_files/seurat_mo8/seurat_8mo_annots_new_0306.rds")
DimPlot(astro_mo8, reduction = "umap", label = T, group.by = "cell_type")
FeaturePlot(astro_mo8, features = "GAD1", reduction = "umap")

# GABAergic inhibitory neurons
INH <- subset(astro_mo8, subset = cell_type == "GABA_INH")
VlnPlot(INH, features = c("KCNN1", "KCNN2"), group.by = "diagnosis") #initial check using violin plots

# GAD1+ cells
gad1_expr <- FetchData(astro_mo8, vars = "GAD1")
ggplot(gad1_expr, aes(x = GAD1)) +
  geom_density(fill = "blue", alpha = 0.5) +
  theme_minimal() +
  labs(title = "GAD1 Expression Density",
       x = "GAD1 Expression", y = "Density")
gad1_cells <- subset(astro_mo8, subset = GAD1 > 0.2) # 0.2 as threshold

# GAD2+ cells
gad2_expr <- FetchData(astro_mo8, vars = "GAD2")
ggplot(gad2_expr, aes(x = GAD2)) +
  geom_density(fill = "pink", alpha = 0.5) +
  theme_minimal() +
  labs(title = "GAD2 Expression Density",
       x = "GAD2 Expression", y = "Density") 
gad2_cells <- subset(astro_mo8, subset = GAD2 > 0.3) # 0.3 as threshold

# Function to create DESeq2 dds object
prepare_dds <- function(seurat) {
  p <- seurat@meta.data
  ids <- table(p$patient_id) # donor is the column name with the cell-line ID
  ids <- names(ids[ids >= 10]) # Select donors with at least 10 cells
  p <- p[p$patient_id %in% ids, ]
  counts <- seurat@assays$RNA$counts #raw counts
  counts <- counts[, rownames(p)]
  sce <- SingleCellExperiment(assays = list(counts = counts), colData = p[, c("patient_id"), drop = FALSE])
  groups <- SingleCellExperiment::colData(sce)
  pb <- Matrix.utils::aggregate.Matrix(t(counts(sce)), 
        groupings = groups, fun = "sum") # Aggregate across cluster-sample groups
  df <- as.data.frame(t(pb))
  p <- p[!duplicated(p$patient_id), ] # Remove duplicates so that p contains only a unique row for each donor
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
  return(list(df = df, p = p, diagnosis = p$diagnosis))
}

# pre-defined GABAergic Inhibitatory Neuons
INH_result <- prepare_dds(INH)
print(table(INH_result$diagnosis))
table(INH_result$p$diagnosis, INH_result$p$sex)
dds_INH <- DESeqDataSetFromMatrix(INH_result$df, colData = INH_result$p, design = ~ diagnosis + sex)
smallestGroupSize <- 2 
keep_INH <- rowSums(counts(dds_INH) >= 10) >= smallestGroupSize 
dds_INH <- dds_INH[keep_INH,] # keep genes with at least 10 counts in two samples
dds_INH <- DESeq(dds_INH)
res_INH <- as.data.frame(results(dds_INH, name = "diagnosis_AD_vs_NC"))
res_INH[c("KCNN1", "KCNN2", "KCNN3", "KCNN4"), ]

# GAD1+ cells
gad1_result <- prepare_dds(gad1_cells)
print(table(gad1_result$diagnosis))
table(gad1_result$p$diagnosis, gad1_result$p$sex)
dds_gad1 <- DESeqDataSetFromMatrix(gad1_result$df, colData = gad1_result$p, design = ~ diagnosis + sex)
smallestGroupSize <- 2 
keep_gad1 <- rowSums(counts(dds_gad1) >= 10) >= smallestGroupSize 
dds_gad1 <- dds_gad1[keep_gad1,] # keep genes with at least 10 counts in two samples
dds_gad1 <- DESeq(dds_gad1)
res_gad1 <- as.data.frame(results(dds_gad1, name = "diagnosis_AD_vs_NC"))
res_gad1[c("KCNN1", "KCNN2", "KCNN3", "KCNN4"), ]

# GAD2+ cells
gad2_result <- prepare_dds_no_celltype(gad2_cells)
print(table(gad2_result$diagnosis))
table(gad2_result$p$diagnosis, gad2_result$p$sex)
dds_gad2 <- DESeqDataSetFromMatrix(gad2_result$df,colData = gad2_result$p, design = ~ diagnosis + sex)
smallestGroupSize <- 2 
keep_gad2 <- rowSums(counts(dds_gad2) >= 10) >= smallestGroupSize 
dds_gad2 <- dds_gad2[keep_gad2,] # keep genes with at least 10 counts in two samples
dds_gad2 <- DESeq(dds_gad2)
res_gad2 <- as.data.frame(results(dds_gad2, name = "diagnosis_AD_vs_NC"))
res_gad2[c("KCNN1", "KCNN2", "KCNN3", "KCNN4"), ]

# Function to plot comparison results
plot_gene_expression <- function(genes_of_interest, dds, title = "", subtitle = "") {
  norm_counts <- as.data.frame(counts(dds, normalized = TRUE)) %>%
    rownames_to_column("gene") %>%
    filter(gene %in% genes_of_interest) %>%
    pivot_longer(cols = -gene, names_to = "sample", values_to = "expression")  
  # Convert colData from DESeq2 object to a metadata dataframe
  meta_df <- as.data.frame(colData(dds)) %>%
    rownames_to_column("sample")
  # Join normalized counts with metadata
  joined_df <- left_join(norm_counts, meta_df, by = "sample")
  # Filter and clean up dataframe for plotting
  plot_df <- joined_df %>%
    dplyr::filter(gene %in% genes_of_interest) %>%
    dplyr::select(all_of(c("gene", "diagnosis", "expression", "patient_id"))) %>%
    dplyr::rename(cell_line = patient_id)
  # Plot
  p <- ggplot(plot_df, aes(x = diagnosis, y = expression, color = cell_line, group = cell_line)) +
    geom_point(size = 3) +
    facet_wrap(~ gene, nrow = 2) +
    theme_bw() +
    labs(title = title, subtitle = subtitle, y = "Normalized RNA-seq expression", x = "Diagnosis", color = "Cell Line") +
    theme(
      text = element_text(size = 14),
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold", size = 16),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, face = "italic", hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.key.size = unit(0.7, "cm"))
  return(p)
}

genes_of_interest <- c("KCNN1", "KCNN2", "KCNN3", "KCNN4")
plot_gene_expression(genes_of_interest, dds_INH, title = "SK channels expression in GABAergic inhibitory neurons",
  subtitle = "8-month organoids")
plot_gene_expression(genes_of_interest, dds_gad1, title = "SK channels expression in GAD1+ Cells",
  subtitle = "8-month organoids")
plot_gene_expression(genes_of_interest, dds_gad2, title = "SK channels expression in GAD2+ cells",
  subtitle = "8-month organoids")
