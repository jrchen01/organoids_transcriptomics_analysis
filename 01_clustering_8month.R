# QC, Harmony integration, clustering, annotation, and Figure 1 - 8-month organoids

library(Seurat)
library(dplyr)
library(data.table)
library(harmony)
library(ggplot2)
library(patchwork)
library(circlize)
library(ComplexHeatmap)

raw_data_path <- "/netapp/LOG-G4/mcuoco/sarah_organoid/scrnaseq_cellranger/preprocess/01_filter/filtered.rds"

## QC filtering, 8-month channel subset
seurat <- readRDS(raw_data_path)

a <- seurat@meta.data
a <- a[a$nFeature_RNA > 200 & a$nFeature_RNA < 6000 & a$percent_mito < 20 & a$`Droplet Type` %in% "SNG", ]
a <- a[a$Channel %in% c("20240321_03112024_4", "20250609_SF-1", "20251001_SF-1"), ]  # 8-month, 19,282 cells

seurat <- subset(seurat, cells = rownames(a))

## normalization, per-channel HVG union, ribosomal module score for regression
seurat <- seurat %>% NormalizeData(verbose = FALSE)

VariableFeatures(seurat) <- split(row.names(seurat@meta.data), seurat@meta.data$Channel) %>%
  lapply(function(cells_use) {
    seurat[, cells_use] %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
      VariableFeatures()
  }) %>% unlist() %>% unique()

translation_genes <- c(
  "FAU", "EEF1A1", "NACA", "RACK1", "RPL10", "RPL10A", "RPL11", "RPL12",
  "RPL13", "RPL13A", "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", "RPL19",
  "RPL21", "RPL22", "RPL23", "RPL23A", "RPL24", "RPL26", "RPL27", "RPL27A",
  "RPL28", "RPL29", "RPL3", "RPL30", "RPL31", "RPL32", "RPL34", "RPL35",
  "RPL35A", "RPL36", "RPL36A", "RPL37", "RPL37A", "RPL38", "RPL39", "RPL4",
  "RPL41", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL8", "RPL9", "RPS10", "RPS11",
  "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", "RPS18", "RPS19",
  "RPS2", "RPS21", "RPS23", "RPS24", "RPS25", "RPS27", "RPS27A", "RPS28",
  "RPS29", "RPS3", "RPS3A", "RPS4X", "RPS5", "RPS6", "RPS7", "RPS8", "RPS9",
  "RPSA", "RPLP0", "RPLP1", "RPLP2", "UBA52"
)
seurat <- AddModuleScore(seurat, features = list(translation_genes), name = "GOBP_Translation")

## scale (regress out translation module score), PCA, Harmony, UMAP
seurat <- seurat %>%
  ScaleData(verbose = FALSE, vars.to.regress = "GOBP_Translation1") %>%
  RunPCA(features = VariableFeatures(seurat), npcs = 30, verbose = FALSE)

seurat <- RunHarmony(seurat, c("Channel", "patient_id"))

seurat <- FindNeighbors(seurat, dims = 1:30, reduction = "harmony", graph.name = "neighbors.harmony.by.donor")
seurat <- RunUMAP(seurat, dims = 1:30, reduction = "harmony", reduction.name = "umap")

## major cell-class clustering and annotation
seurat <- FindClusters(seurat, resolution = 0.05, graph.name = "neighbors.harmony.by.donor")

s.genes <- fread("~/data1/SarahF/organoids/G1.S.cellCycle.genes", data.table = FALSE, header = FALSE)$V1
g2m.genes <- fread("~/data1/SarahF/organoids/G2.M.cellCycle.genes", data.table = FALSE, header = FALSE)$V1
seurat <- CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

seurat$cell <- factor(seurat$seurat_clusters)
levels(seurat$cell) <- c("Glial cells", "Inhibitory neurons", "Excitatory neurons",
                          "Dividing cells", "Fibroblast", "Excitatory neurons")
Idents(seurat) <- "cell"

## subcluster-level clustering and annotation
# moved the diagnosis merge below to right after this block - in the original
# interactive session it happened much later, but seurat@meta.data gets fully
# overwritten there and needs to keep the subcluster column
# cluster "0" here is from the resolution=0.3 run right above it, taking that on faith from the source
seurat <- FindClusters(seurat, resolution = 0.3, graph.name = "neighbors.harmony.by.donor")

subcluster_result <- FindSubCluster(seurat, cluster = "0",
                                     graph.name = "neighbors.harmony.by.donor",
                                     subcluster.name = "subcluster",
                                     resolution = 0.2)

subcluster_result$subcluster <- factor(subcluster_result$subcluster)
levels(subcluster_result$subcluster) <- c("Inhibitory 1", "Inhibitory 2", "Mixed cells", "Dividing cells",
                                           "Excitatory 1", "Astrocyte", "OPC 1", "OPC 2", "Fibroblast",
                                           "Excitatory 2", "Unknown")
seurat$subcluster <- subcluster_result$subcluster
seurat$subcluster <- factor(seurat$subcluster, levels = c("Inhibitory 1", "Inhibitory 2", "Excitatory 1",
                                                            "Excitatory 2", "Dividing cells", "Mixed cells",
                                                            "OPC 1", "OPC 2", "Astrocyte", "Fibroblast", "Unknown"))
Idents(seurat) <- "subcluster"

## diagnosis assignment
# donor IDs show up in a couple different spellings (UCI22 vs UCI-22 etc), matching both
patient <- seurat@meta.data
patient <- patient[!duplicated(patient$patient_id), ]
patient$diagnosis <- NA
patient[patient$patient_id %in% c("UCI22", "UCI-22", "ADRC_2800", "ADRC-2800", "UCI2", "UCI-2",
                                   "UCI96", "UCI-96", "2991", "ADRC-2991", "3341", "ADRC-3341"), "diagnosis"] <- "sAD"
patient[patient$patient_id %in% c("UCI52", "UCI-52", "UCI40", "UCI-40c7", "3551", "ADRC-3551",
                                   "40", "ADRC-40", "3483", "ADRC-3483"), "diagnosis"] <- "CTRL"

p <- seurat@meta.data
p$cellID <- rownames(p)
p <- merge(p, patient[, c("patient_id", "diagnosis")], by = "patient_id")
rownames(p) <- p$cellID
p <- p[colnames(seurat), ]
seurat@meta.data <- p

save(seurat, file = "/raidixshare_logg01/thais/SarahF/organoids/seurat_8mo_harmony_regressedRibog.rda")

# Figure 1: UMAP + marker heatmaps
# meta.data$cell = major class, meta.data$subcluster = subtype

df_toPlot <- cbind(Embeddings(object = seurat, reduction = "umap"), seurat@meta.data)
colnames(df_toPlot)[c(1:2)] <- c("UMAP1", "UMAP2")

## color palettes
cols_main <- c(
  "Glial cells"        = "#d8b7dd",
  "Inhibitory neurons" = "#4f81bd",
  "Excitatory neurons" = "#6d8c3c",
  "Dividing cells"     = "#5abf4b",
  "Fibroblast"         = "#b2182b"
)

cols_sub <- c(
  "#1f4e79",  # Inhibitory 1
  "#4f81bd",  # Inhibitory 2
  "#6d8c3c",  # Excitatory 1
  "#3b6b2f",  # Excitatory 2
  "#5abf4b",  # Dividing
  "#e6c8a8",  # Mixed
  "#8c510a",  # OPC 1
  "#d8b7dd",  # OPC 2
  "#7b3294",  # Astrocyte
  "#b2182b",  # Fibroblast
  "grey70"    # Unknown
)

## Fig 1A/1C - UMAP, major class and subcluster
p_sub <- ggplot(df_toPlot, aes(x = UMAP1, y = UMAP2, color = subcluster)) +
  geom_point(size = 0.6, stroke = 0) +
  scale_color_manual(values = cols_sub) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 10),
    legend.key.size = unit(0.8, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))

p_main <- ggplot(df_toPlot, aes(x = UMAP1, y = UMAP2, color = cell)) +
  geom_point(size = 0.6, stroke = 0) +
  scale_color_manual(values = cols_main) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 10),
    legend.key.size = unit(0.8, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))

pdf("/raidixshare_logg01/jchen/project/01_sarah_proj/05_organoids/02_results/05_clustering/clustering_umap/8mo_umap.pdf", width = 12, height = 6)
p_sub + p_main
dev.off()

## Fig 1B/1D - marker heatmaps (z-score of mean expression)
marker_genes <- c("TUBB3", "MAP2", "GAD1", "GAD2", "SLC17A7", "SLC17A6",
                   "SOX9", "CENPK", "OLIG1", "PDGFRA", "GFAP", "AQP4")

col_fun <- circlize::colorRamp2(c(-1, 0, 4), c("white", "white", "red"))

make_marker_heatmap <- function(seurat_obj, group_by, genes) {
  a <- AverageExpression(object = seurat_obj, group.by = group_by)$RNA
  a <- as.data.frame(a)
  aux <- a[genes, colnames(a)]
  scaled_mat <- t(scale(t(aux)))

  Heatmap(scaled_mat,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          col = col_fun,
          name = "Z-score",
          row_names_side = "left",
          row_names_gp = gpar(fontsize = 9, fontface = "italic"),
          column_names_side = "top",
          column_names_rot = 45,
          column_names_gp = gpar(fontsize = 9),
          rect_gp = gpar(col = "grey90", lwd = 0.5),
          border = TRUE,
          border_gp = gpar(col = "black", lwd = 1),
          heatmap_legend_param = list(
            title = "Z-score",
            title_gp = gpar(fontsize = 9, fontface = "bold"),
            labels_gp = gpar(fontsize = 8),
            legend_height = unit(3, "cm"),
            legend_width = unit(0.4, "cm")
          ),
          show_row_dend = FALSE,
          show_column_dend = FALSE)
}

heatmap_8mo_sub  <- make_marker_heatmap(seurat, "subcluster", marker_genes)
heatmap_8mo_main <- make_marker_heatmap(seurat, "cell", marker_genes)

pdf("/raidixshare_logg01/jchen/project/01_sarah_proj/05_organoids/02_results/05_clustering/clustering_umap/8mo_heatmap.pdf", width = 8, height = 8)
grid::grid.newpage()
draw(heatmap_8mo_main, padding = unit(c(10, 10, 5, 10), "mm"), newpage = FALSE)
grid::grid.newpage()
draw(heatmap_8mo_sub, padding = unit(c(10, 10, 5, 10), "mm"), newpage = FALSE)
dev.off()
