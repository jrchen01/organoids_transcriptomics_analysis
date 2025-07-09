#!/usr/bin/env Rscript

# Cell type annotation
#devtools::install_github("immunogenomics/presto")
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(devtools)
library(presto)

# Step1 - markers for each cluster 
## Step1.1 - check *non-unique markers* for each cluster 
astro_mo8 <- readRDS("/raidixshare_logg01/jchen/project/01_sarah_proj/03_seurat_files/seurat_mo8/seurat_8mo_dim30_res0.15.rds")
DimPlot(astro_mo8, reduction = "umap") #check
head(astro_mo8@meta.data, 3) #check
Idents(astro_mo8) <- "RNA_snn_res.0.15"
# Function to plot a volcano plot using negLOG10P and LOG2FC
plot_volcano <- function(df) {
  # Convert necessary columns to numeric
  df$LOG2FC <- as.numeric(df$avg_log2FC)  # Use avg_log2FC as log2 fold change
  df$negLOG10P <- as.numeric(-log10(df$p_val_adj))  # Use p_val_adj for negLOG10P
  # Handle missing row names (gene names)
  df$gene <- rownames(df)  # Create a gene column from rownames
  # Categorize genes based on log2 fold change and p-value
  df <- df %>%
    mutate(type = case_when(
      LOG2FC > 2 & negLOG10P > 3 ~ "P < 0.001, log2FC > 2", 
      LOG2FC < -2 & negLOG10P > 3 ~ "P < 0.001, log2FC < -2", 
      TRUE ~ "None"
    )) %>%
    mutate(type = factor(type, levels = c("P < 0.001, log2FC < -2", "None", "P < 0.001, log2FC > 2")))
  # Create volcano plot
  ggplot(df, aes(LOG2FC, negLOG10P)) +
    geom_point(aes(color = type), size = 1.5) + 
    scale_color_manual(values = c("blue", "gray", "red"), name = NULL) + 
    geom_hline(yintercept = -log10(0.05), linetype = 2) + 
    geom_hline(yintercept = -log10(0.001), linetype = 2) + 
    geom_vline(xintercept = c(-2, 2), linetype = 2) + 
    geom_text_repel(
      data = df[df$negLOG10P > 3 & abs(df$LOG2FC) > 2,],
      aes(label = gene),  # Use the new gene column for labeling
      size = 2, color = "black", segment.color = "black",
      fontface = "bold", max.overlaps = 20, show.legend = FALSE
    ) +
    theme_bw() +
    theme(legend.position.inside = c(0.85, 0.15),  # Use the new legend position argument
          legend.background = element_blank())
}
# Example usage:
# plot_volcano(deg_cluster2)

for (i in 0:7) {
  assign(paste0("deg_cluster", i), FindMarkers(astro_mo8, ident.1 = i))}
for (i in 0:7) {
  cat("deg_cluster", i, ": ", nrow(get(paste0("deg_cluster", i))), "\n")}

plot_volcano(deg_cluster0)
deg_cluster0_down <- deg_cluster0 %>%
  filter(avg_log2FC < -2, p_val_adj < 0.001)
nrow(deg_cluster0_down) #470
deg_cluster0_up <- deg_cluster0 %>%
  filter(avg_log2FC > 1, p_val_adj < 0.001)
nrow(deg_cluster0_up) #60
genes_cluster0_down <- rownames(deg_cluster0_down)
write.table(genes_cluster0_down, file = "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/plots/month8_seurat/01_dim30_res0.15/DEGs/deg_cluster0_down_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
genes_cluster0_up <- rownames(deg_cluster0_up)
write.table(genes_cluster0_up, file = "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/plots/month8_seurat/01_dim30_res0.15/DEGs/deg_cluster0_up_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

plot_volcano(deg_cluster1)
deg_cluster1_down <- deg_cluster1 %>%
  filter(avg_log2FC < -1.5, p_val_adj < 0.001)
nrow(deg_cluster1_down) #148
deg_cluster1_up <- deg_cluster1 %>%
  filter(avg_log2FC > 1.5, p_val_adj < 0.001)
nrow(deg_cluster1_up) #241
genes_cluster1_down <- rownames(deg_cluster1_down)
write.table(genes_cluster1_down, file = "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/plots/month8_seurat/01_dim30_res0.15/DEGs/deg_cluster1_down_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
genes_cluster1_up <- rownames(deg_cluster1_up)
write.table(genes_cluster1_up, file = "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/plots/month8_seurat/01_dim30_res0.15/DEGs/deg_cluster1_up_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

plot_volcano(deg_cluster2)
deg_cluster2_down <- deg_cluster2 %>%
  filter(avg_log2FC < -1.5, p_val_adj < 0.001)
nrow(deg_cluster2_down) #172
deg_cluster2_up <- deg_cluster2 %>%
  filter(avg_log2FC > 1, p_val_adj < 0.001)
nrow(deg_cluster2_up) #69
genes_cluster2_down <- rownames(deg_cluster2_down)
write.table(genes_cluster2_down, file = "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/plots/month8_seurat/01_dim30_res0.15/DEGs/deg_cluster2_down_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
genes_cluster2_up <- rownames(deg_cluster2_up)
write.table(genes_cluster2_up, file = "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/plots/month8_seurat/01_dim30_res0.15/DEGs/deg_cluster2_up_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

plot_volcano(deg_cluster3)
deg_cluster3_down <- deg_cluster3 %>%
  filter(avg_log2FC < -1.5, p_val_adj < 0.001)
nrow(deg_cluster3_down) #85
deg_cluster3_up <- deg_cluster3 %>%
  filter(avg_log2FC > 1.5, p_val_adj < 0.001)
nrow(deg_cluster3_up) #221
genes_cluster3_down <- rownames(deg_cluster3_down)
write.table(genes_cluster3_down, file = "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/plots/month8_seurat/01_dim30_res0.15/DEGs/deg_cluster3_down_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
genes_cluster3_up <- rownames(deg_cluster3_up)
write.table(genes_cluster3_up, file = "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/plots/month8_seurat/01_dim30_res0.15/DEGs/deg_cluster3_up_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

plot_volcano(deg_cluster4)
deg_cluster4_down <- deg_cluster4 %>%
  filter(avg_log2FC < -2, p_val_adj < 0.001)
nrow(deg_cluster4_down) #167
deg_cluster4_up <- deg_cluster4 %>%
  filter(avg_log2FC > 2, p_val_adj < 0.001)
nrow(deg_cluster4_up) #339
genes_cluster4_down <- rownames(deg_cluster4_down)
write.table(genes_cluster4_down, file = "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/plots/month8_seurat/01_dim30_res0.15/DEGs/deg_cluster4_down_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
genes_cluster4_up <- rownames(deg_cluster4_up)
write.table(genes_cluster4_up, file = "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/plots/month8_seurat/01_dim30_res0.15/DEGs/deg_cluster4_up_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

plot_volcano(deg_cluster5)
deg_cluster5_down <- deg_cluster5 %>%
  filter(avg_log2FC < -2, p_val_adj < 0.001)
nrow(deg_cluster5_down) #439
deg_cluster5_up <- deg_cluster5 %>%
  filter(avg_log2FC > 2, p_val_adj < 0.001)
nrow(deg_cluster5_up) #868
genes_cluster5_down <- rownames(deg_cluster5_down)
write.table(genes_cluster5_down, file = "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/plots/month8_seurat/01_dim30_res0.15/DEGs/deg_cluster5_down_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
genes_cluster5_up <- rownames(deg_cluster5_up)
write.table(genes_cluster5_up, file = "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/plots/month8_seurat/01_dim30_res0.15/DEGs/deg_cluster5_up_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

plot_volcano(deg_cluster6)
deg_cluster6_down <- deg_cluster6 %>%
  filter(avg_log2FC < -2, p_val_adj < 0.001)
nrow(deg_cluster6_down) #265
deg_cluster6_up <- deg_cluster6 %>%
  filter(avg_log2FC > 2, p_val_adj < 0.001)
nrow(deg_cluster6_up) #771
genes_cluster6_down <- rownames(deg_cluster6_down)
write.table(genes_cluster6_down, file = "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/plots/month8_seurat/01_dim30_res0.15/DEGs/deg_cluster6_down_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
genes_cluster6_up <- rownames(deg_cluster6_up)
write.table(genes_cluster6_up, file = "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/plots/month8_seurat/01_dim30_res0.15/DEGs/deg_cluster6_up_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

plot_volcano(deg_cluster7)
deg_cluster7_down <- deg_cluster7 %>%
  filter(avg_log2FC < -2, p_val_adj < 0.001)
nrow(deg_cluster7_down) #286
deg_cluster7_up <- deg_cluster7 %>%
  filter(avg_log2FC > 2.5, p_val_adj < 1e-100)
nrow(deg_cluster7_up) #341
genes_cluster7_down <- rownames(deg_cluster7_down)
write.table(genes_cluster7_down, file = "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/plots/month8_seurat/01_dim30_res0.15/DEGs/deg_cluster7_down_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
genes_cluster7_up <- rownames(deg_cluster7_up)
write.table(genes_cluster7_up, file = "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/plots/month8_seurat/01_dim30_res0.15/DEGs/deg_cluster7_up_genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

## Step1.2 - check *unique markers* for each cluster
### save all non-unique markers
base_path <- "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/Seurat_8mo/01_degs/01_markers_each_cluster/"
for (i in 0:7) {
  deg_data <- get(paste0("deg_cluster", i))
  file_path <- paste0(base_path, "cluster_", i, "_markers.csv")
  write.csv(deg_data, file = file_path, row.names = FALSE)}

for (i in 0:7) {
  deg_data <- get(paste0("deg_cluster", i))
  flt_deg_data <- deg_data %>%
    filter(abs(avg_log2FC) > 1, p_val_adj < 0.001) # filtering cutoff: log2FC=1; p=0.001
  assign(paste0("flt_deg_cluster", i), flt_deg_data)}
rows_list <- lapply(0:7, function(i) {
  flt_deg_data <- get(paste0("flt_deg_cluster", i))
  nrow(flt_deg_data)})
names(rows_list) <- paste0("flt_deg_cluster", 0:7)
print(rows_list)

for (i in 0:7) {
  deg_data <- get(paste0("flt_deg_cluster", i))
  deg_data$gene <- rownames(deg_data)
  assign(paste0("flt_deg_cluster", i), deg_data)}

deg_clusters <- list(
  deg_cluster0 = flt_deg_cluster0,
  deg_cluster1 = flt_deg_cluster1,
  deg_cluster2 = flt_deg_cluster2,
  deg_cluster3 = flt_deg_cluster3,
  deg_cluster4 = flt_deg_cluster4,
  deg_cluster5 = flt_deg_cluster5,
  deg_cluster6 = flt_deg_cluster6,
  deg_cluster7 = flt_deg_cluster7)

unique_markers <- lapply(names(deg_clusters), function(cluster_name) {
  current_deg <- deg_clusters[[cluster_name]]
  current_genes <- current_deg$gene
  other_genes <- unlist(lapply(deg_clusters[names(deg_clusters) != cluster_name], function(df) df$gene))
  unique_genes <- setdiff(current_genes, other_genes)
  current_unique <- current_deg[current_deg$gene %in% unique_genes, ]
  return(current_unique)})
names(unique_markers) <- paste0("unique_deg_cl", 0:7)

for (i in 0:7) {
  name <- paste0("unique_deg_cl", i)
  assign(name, unique_markers[[name]])}
                        
#### cluster 0
nrow(unique_deg_cl0) #84
unique_deg_cluster0_down <- unique_deg_cl0 %>%
  filter(avg_log2FC < -1)
nrow(unique_deg_cluster0_down) #59
unique_deg_cluster0_up <- unique_deg_cl0 %>%
  filter(avg_log2FC > 1)
nrow(unique_deg_cluster0_up) #25
#### cluster 1
nrow(unique_deg_cl1) #176
unique_deg_cluster1_down <- unique_deg_cl1 %>%
  filter(avg_log2FC < -1)
nrow(unique_deg_cluster1_down) #12
unique_deg_cluster1_up <- unique_deg_cl1 %>%
  filter(avg_log2FC > 1)
nrow(unique_deg_cluster1_up) #164
#### Cluster 2
nrow(unique_deg_cl2) #365
unique_deg_cluster2_down <- unique_deg_cl2 %>%
  filter(avg_log2FC < -1)
nrow(unique_deg_cluster2_down) #356
unique_deg_cluster2_up <- unique_deg_cl2 %>%
  filter(avg_log2FC > 1)
nrow(unique_deg_cluster2_up) #9
unique_deg_cluster2_down <- unique_deg_cluster2_down %>%
  arrange(p_val_adj) %>%
  slice_head(n = 200)
nrow(unique_deg_cluster2_down) #extract the top 200
#### Cluster 3
nrow(unique_deg_cl3) #116
unique_deg_cluster3_down <- unique_deg_cl3 %>%
  filter(avg_log2FC < -1)
nrow(unique_deg_cluster3_down) #2
unique_deg_cluster3_up <- unique_deg_cl3 %>%
  filter(avg_log2FC > 1)
nrow(unique_deg_cluster3_up) #114
#### Cluster 4
nrow(unique_deg_cl4) #308
unique_deg_cluster4_down <- unique_deg_cl4 %>%
  filter(avg_log2FC < -1)
nrow(unique_deg_cluster4_down) #15
unique_deg_cluster4_up <- unique_deg_cl4 %>%
  filter(avg_log2FC > 1)
nrow(unique_deg_cluster4_up) #293
unique_deg_cluster4_up <- unique_deg_cluster4_up %>%
  arrange(p_val_adj) %>%
  slice_head(n = 200)
nrow(unique_deg_cluster4_up) #extract the top 200
#### Cluster 5
nrow(unique_deg_cl5) #917
unique_deg_cluster5_down <- unique_deg_cl5 %>%
  filter(avg_log2FC < -1)
nrow(unique_deg_cluster5_down) #266
unique_deg_cluster5_up <- unique_deg_cl5 %>%
  filter(avg_log2FC > 1)
nrow(unique_deg_cluster5_up) #651
unique_deg_cluster5_down <- unique_deg_cluster5_down %>%
  arrange(p_val_adj) %>%
  slice_head(n = 200)
unique_deg_cluster5_up <- unique_deg_cluster5_up %>%
  arrange(p_val_adj) %>%
  slice_head(n = 200)
nrow(unique_deg_cluster5_down) #extract the top 200
nrow(unique_deg_cluster5_up) #extract the top 200
#### Cluster 6
nrow(unique_deg_cl6) #2032
unique_deg_cluster6_down <- unique_deg_cl6 %>%
  filter(avg_log2FC < -1)
nrow(unique_deg_cluster6_down) #1073
unique_deg_cluster6_up <- unique_deg_cl6 %>%
  filter(avg_log2FC > 1)
nrow(unique_deg_cluster6_up) #959
unique_deg_cluster6_down <- unique_deg_cluster6_down %>%
  arrange(p_val_adj) %>%
  slice_head(n = 200)
unique_deg_cluster6_up <- unique_deg_cluster6_up %>%
  arrange(p_val_adj) %>%
  slice_head(n = 200)
nrow(unique_deg_cluster6_down) #extract the top 200
nrow(unique_deg_cluster6_up) #extract the top 200
#### Cluster 7
nrow(unique_deg_cl7) #1537
unique_deg_cluster7_down <- unique_deg_cl7 %>%
  filter(avg_log2FC < -1)
nrow(unique_deg_cluster7_down) #44
unique_deg_cluster7_up <- unique_deg_cl7 %>%
  filter(avg_log2FC > 1)
nrow(unique_deg_cluster7_up) #1493
unique_deg_cluster7_up <- unique_deg_cluster7_up %>%
  arrange(p_val_adj) %>%
  slice_head(n = 200)
nrow(unique_deg_cluster7_up) #extract the top 200

base_path <- "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/plots/month8_seurat/01_dim30_res0.15/DEGs/01_unique_markers/"
for (i in 0:7) { #save unique markers
  up_genes <- rownames(get(paste0("unique_deg_cluster", i, "_up")))
  down_genes <- rownames(get(paste0("unique_deg_cluster", i, "_down")))
  up_file <- paste0(base_path, "unique_deg_cluster", i, "_up.txt")
  down_file <- paste0(base_path, "unique_deg_cluster", i, "_down.txt")
  write.table(up_genes, file = up_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(down_genes, file = down_file, row.names = FALSE, col.names = FALSE, quote = FALSE)}

rm(list=ls())

# Step2 - add module score
library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

astro_mo8 <- readRDS("/raidixshare_logg01/jchen/project/01_sarah_proj/03_seurat_files/seurat_mo8/seurat_8mo_dim30_res0.15.rds")
DimPlot(astro_mo8, reduction = "umap") #just check


plot_module_score <- function(seurat_obj, marker_list, module_name, output_path) {
  seurat_obj <- AddModuleScore(
    object = seurat_obj,
    features = marker_list,
    name = module_name)
  
  for (i in 1:length(marker_list)) {
    feature_name <- paste0(module_name, i)
    cell_type <- names(marker_list)[i]
    output_file <- paste0(output_path, "moduleScore_", gsub(" ", "_", tolower(cell_type)), ".png")
    
    p <- FeaturePlot(
      object = seurat_obj,
      features = feature_name,
      cols = rev(brewer.pal(n = 11, name = "RdBu")),
      label = TRUE,
      reduction = "umap"
    ) +
      ggtitle(cell_type) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    
    ggsave(filename = output_file, plot = p, width = 8, height = 6, dpi = 300)
  }
  return(seurat_obj)
}

output_path <- "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/Seurat_8mo/plots/01_moduleScore_plots/updated/" 

marker_genes <- list(
  Cell_Death = c("CASP3", "CASP9", "BAX", "BCL2", "TP53", "CYCS", "FAS", "RIPK1", "RIPK3", "MLKL"),
  Cell_Cycle = c("CDK1", "CCNB1", "CCNA2", "CDK2", "CCND1", "PCNA", "MCM2", "MCM6", "CDC25A", "AURKA", "AURKB"),
  Cell_Stress = c("HMOX1", "NQO1", "SOD2", "HSPA5", "ATF4", "ATF6", "HIF1A", "VEGFA", "LDHA", "TP53", "ATM", "ATR", "IL6", "TNF", "SAA1", "CCL2"),
  Autophagy = c("BECN1", "ATG5", "ATG7"),
  Inhibitory_Neurons = c("GAD1", "GAD2", "SLC6A1", "GABBR1", "GABRA1", "LHX6", "DLX1", "DLX2", "SST", "PVALB", "NPY", "VIP", "CNR1", "RELN", "SNAP25", "ARX", "CRHBP", "COCH", "SYNPR"),
  Excitatory_Neurons = c("SLC17A7", "SLC17A6", "SLC1A1", "GRIN1", "GRIN2B", "SATB2", "NRXN3", "SYT1", "SLC17A7", "NLGN1", "NEUROD6", "TBR1", "BCL11B"),
  OPC_ODC = c("PDGFRA", "OLIG2", "SOX10", "CSPG4", "MYT1", "MBP", "MAG", "MOG", "PLP1", "APOD", "GPR17", "FABP7", "QKI", "BCAS1", "TSPAN2"),
  Astrocyte = c("SLC1A3", "SLC1A2", "GLUL", "SOX9", "AQP4", "ALDH1L1", "GFAP", "GJA1", "S100B", "CD44", "PTGDS", "NDRG2", "VIM", "HOPX", "SERPINA3"),
  Glial_Cell_General = c("GFAP", "S100B", "SLC1A3", "ALDH1L1", "AQP4", "GJA1", "SOX9", "CD44", "CNP", "PLP1", "MBP", "MAG", "OLIG2", "PDGFRA", "TNC", "CLU", "SERPINA3", "SOX10", "CSPG4", "MYT1", "APOD"),
  Mature_Neurons = c("MAP2", "RBFOX3", "SYT1", "TUBB3", "SNAP25", "GRIN1", "GRIN2B", "SATB2", "TBR1", "CAMK2A", "SLC17A7", "SLC17A6", "BCL11B", "FOXP1", "STMN2", "NEFH"),
  Neuron_Progenitors = c("SOX2", "NES", "PAX6", "VIM", "ASCL1", "HES1", "FABP7", "DCX", "EOMES", "NOTCH1", "ID4", "CD24", "SOX4", "SOX11"),
  Immature_Neurons = c("DCX", "EOMES", "ASCL1", "NEUROD1", "NEUROG2", "STMN1", "DCLK1", "SOX4", "SOX11", "PROX1", "DLX2", "CD24", "GAP43", "DLX1", "INSM1", "NHLH1")
)

plot_module_score(astro_mo8, marker_genes, "Updated_Markers", output_path)

rm(list=ls())

# additional analysis on NKCC1 and KCC2 expression (in cluster 5 and 6)
library(Seurat)
library(dplyr)
library(ggplot2)

astro_mo8 <- readRDS("/raidixshare_logg01/jchen/project/01_sarah_proj/03_seurat_files/seurat_mo8/seurat_8mo_dim30_res0.15.rds")
Idents(astro_mo8) <- "seurat_clusters"
candidate_genes <- c("SLC12A2", "SLC12A5")
target_cells <- WhichCells(astro_mo8, idents = "5")
expr_mat <- GetAssayData(astro_mo8, slot = "data", assay = "RNA")
avg_expr <- rowMeans(expr_mat[, target_cells, drop = FALSE]) 
for (g in candidate_genes) { #histogram
  hist(avg_expr, xlim = c(0, 2), breaks = 50, main = paste("cluster_5", "\nGene:", g), xlab = "Average log-normalized expression")
  if (g %in% names(avg_expr)) {
    abline(v = avg_expr[g], col = "red", lty = 2, lwd = 2)
    text(x = avg_expr[g], y = 0.6 * max(hist(avg_expr, plot = FALSE)$counts),
         labels = g, col = "red", srt = 90, pos = 4)}}

target_cells <- WhichCells(astro_mo8, idents = "6")
expr_mat <- GetAssayData(astro_mo8, slot = "data", assay = "RNA")
avg_expr <- rowMeans(expr_mat[, target_cells, drop = FALSE]) 
for (g in candidate_genes) { #histogram
  hist(avg_expr, xlim = c(0, 2), breaks = 50, main = paste("cluster_6", "\nGene:", g), xlab = "Average log-normalized expression")
  if (g %in% names(avg_expr)) {
    abline(v = avg_expr[g], col = "red", lty = 2, lwd = 2)
    text(x = avg_expr[g], y = 0.6 * max(hist(avg_expr, plot = FALSE)$counts),
         labels = g, col = "red", srt = 90, pos = 4)}}

rm(list=ls())
# cell type annotation
library(Seurat)
library(dplyr)
library(ggplot2)

## subcluster and additional analysis
### Subcluster for the glia cell cluster - cluster 4
astro_mo8 <- readRDS("/raidixshare_logg01/jchen/project/01_sarah_proj/03_seurat_files/seurat_mo8/seurat_8mo_dim30_res0.15.rds")
DimPlot(astro_mo8, reduction = "umap")
names(astro_mo8@graphs) # check
astro_mo8 <- FindSubCluster(object = astro_mo8, cluster = 4, graph.name = "RNA_snn",     
  subcluster.name = "sub.cluster", resolution = 0.5, algorithm = 1)

DimPlot(astro_mo8, reduction = "umap", group.by = "sub.cluster", label = TRUE) +
  ggtitle("Subcluster")
Idents(astro_mo8) <- "sub.cluster"
cells_4_0 <- WhichCells(astro_mo8, idents = "4_0")
astro_mo8@meta.data[cells_4_0, "cell_type"] <- "OPC"
cells_4_1 <- WhichCells(astro_mo8, idents = "4_1")
astro_mo8@meta.data[cells_4_1, "cell_type"] <- "ASC"
cells_4_2 <- WhichCells(astro_mo8, idents = "4_2")
astro_mo8@meta.data[cells_4_2, "cell_type"] <- "TBD"
cells_4_3 <- WhichCells(astro_mo8, idents = "4_3")
astro_mo8@meta.data[cells_4_3, "cell_type"] <- "TBD"

Idents(astro_mo8) <- "seurat_clusters"
cells_0 <- WhichCells(astro_mo8, idents = "0")
astro_mo8@meta.data[cells_0, "cell_type"] <- "TBD" # might be some neurons
cells_1 <- WhichCells(astro_mo8, idents = "1")
astro_mo8@meta.data[cells_1, "cell_type"] <- "Dividing_cells"
cells_2 <- WhichCells(astro_mo8, idents = "2")
astro_mo8@meta.data[cells_2, "cell_type"] <- "TBD"
cells_3 <- WhichCells(astro_mo8, idents = "3")
astro_mo8@meta.data[cells_3, "cell_type"] <- "Dividing_cells"
cells_5 <- WhichCells(astro_mo8, idents = "5")
astro_mo8@meta.data[cells_5, "cell_type"] <- "EX"
cells_6 <- WhichCells(astro_mo8, idents = "6")
astro_mo8@meta.data[cells_6, "cell_type"] <- "NKCC1_EX"
cells_7 <- WhichCells(astro_mo8, idents = "7")
astro_mo8@meta.data[cells_7, "cell_type"] <- "TBD"
DimPlot(astro_mo8, group.by = "cell_type", reduction = "umap", label = T, repel = TRUE) #just check

### Subcluster for the neuron cluster - cluster 0
astro_mo8 <- FindSubCluster(object = astro_mo8, cluster = 0, graph.name = "RNA_snn",     
  subcluster.name = "sub.cluster", resolution = 0.5, algorithm = 1)
FeaturePlot(astro_mo8, features = "GAD1", reduction = "umap", cols = c("lightgrey", "blue")) # GAD1 expression
FeaturePlot(astro_mo8, features = "SLC6A1", reduction = "umap", cols = c("lightgrey", "red")) #GAT1 expression
Idents(astro_mo8) <- "sub.cluster"
cells_0_0 <- WhichCells(astro_mo8, idents = "0_0")
astro_mo8@meta.data[cells_0_0, "cell_type"] <- "INH"
DimPlot(astro_mo8, group.by = "cell_type", reduction = "umap", label = T, repel = TRUE) #check
table(astro_mo8@meta.data$cell_type) #check

saveRDS(astro_mo8, "/raidixshare_logg01/jchen/project/01_sarah_proj/03_seurat_files/seurat_mo8/seurat_8mo_annots_0411.rds")

# Additional analysis in "GABA_INH" cluster - check GABAergic marker expression
rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
astro_mo8 <- readRDS("/raidixshare_logg01/jchen/project/01_sarah_proj/03_seurat_files/seurat_mo8/seurat_8mo_annots_0411.rds")
Idents(astro_mo8) <- "cell_type"
candidate_genes <- c("GAD1", "GAD2", "SLC6A1", "SLC32A1") # GAD1, GAD2, VGAT (SLC32A1), GAT‑1 (SLC6A1)
target_cells <- WhichCells(astro_mo8, idents = "GABA_INH")
expr_mat <- GetAssayData(astro_mo8, slot = "data", assay = "RNA")
avg_expr <- rowMeans(expr_mat[, target_cells, drop = FALSE]) 
for (g in candidate_genes) { # histogram
  hist(avg_expr, xlim = c(0, 2), breaks = 50, main = paste("GABA_INH", "\nGene:", g), xlab = "Average log-normalized expression")
  if (g %in% names(avg_expr)) {
    abline(v = avg_expr[g], col = "red", lty = 2, lwd = 2)
    text(x = avg_expr[g], y = 0.6 * max(hist(avg_expr, plot = FALSE)$counts),
         labels = g, col = "red", srt = 90, pos = 4)}}
