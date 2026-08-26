# QC, Harmony integration, clustering, annotation, and Figure 2 - 11-month (>10 to 13 months) organoids
## meta.data$cell = subcluster label at this timepoint, meta.data$cell.major = major class (flipped naming vs 8-month object)

library(Seurat)
library(dplyr)
library(data.table)
library(harmony)
library(ggplot2)
library(patchwork)
library(circlize)
library(ComplexHeatmap)

raw_data_path <- "/netapp/LOG-G4/mcuoco/sarah_organoid/scrnaseq_cellranger/preprocess/01_filter/filtered.rds"

## QC filtering, 11-month channel subset
# 3483 and ADRC_2800 were part of the 8-month cohort but aren't in the 11-month one
seurat <- readRDS(raw_data_path)

a <- seurat@meta.data
a <- a[a$nFeature_RNA > 200 & a$nFeature_RNA < 6000 & a$percent_mito < 20 & a$`Droplet Type` %in% "SNG", ]
a <- a[a$Channel %in% c("20250316_SF-1", "20250316_SF-2", "20250316_SF-3", "20251001_SF-2"), ]
a <- a[!a$patient_id %in% c("3483", "ADRC_2800"), ]

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

## major cell-class clustering and annotation (17 clusters at res 0.4)
seurat <- FindClusters(seurat, resolution = 0.4, graph.name = "neighbors.harmony.by.donor")

s.genes <- fread("~/data1/SarahF/organoids/G1.S.cellCycle.genes", data.table = FALSE, header = FALSE)$V1
g2m.genes <- fread("~/data1/SarahF/organoids/G2.M.cellCycle.genes", data.table = FALSE, header = FALSE)$V1
seurat <- CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

seurat$cell <- factor(seurat$seurat_clusters)
seurat$cell <- factor(seurat$cell, levels = 0:16)
levels(seurat$cell) <- c("Astrocyte 1", "Astrocyte 2", "Inhibitory 1", "Dividing NPC", "Fibroblast 1",
                          "Inhibitory 2", "Inhibitory 4", "OPC", "Astrocyte 3", "Excitatory 1",
                          "Senescent cells", "Unknown 1", "Excitatory 2", "Fibroblast 2",
                          "Dedifferentiating cells", "Unknown 2", "Inhibitory 1")
# clusters 2 and 16 both get the "Inhibitory 1" label and merge into one level
Idents(seurat) <- "cell"

## split "Inhibitory 1" into Inhibitory 1 / Inhibitory 3, reclassify cycling
## "Inhibitory 4" cells as dividing cells
subcluster_result <- FindSubCluster(seurat, cluster = "Inhibitory 1",
                                     graph.name = "neighbors.harmony.by.donor",
                                     subcluster.name = "subcluster",
                                     resolution = 0.1)

subcluster_result$subcluster <- factor(subcluster_result$subcluster)
levels(subcluster_result$subcluster)[10:12] <- c("Inhibitory 1", "Inhibitory 3", "Inhibitory 1")
seurat$cell <- subcluster_result$subcluster

p <- seurat@meta.data
p$cell <- as.character(p$cell)
p[p$cell %in% "Inhibitory 4" & !p$Phase %in% "G1", "cell"] <- "Dividing cells"
p$cell <- factor(p$cell, levels = c("Inhibitory 1", "Inhibitory 2", "Inhibitory 3", "Inhibitory 4",
                                     "Excitatory 1", "Excitatory 2", "Dividing NPC", "Dividing cells",
                                     "OPC", "Astrocyte 1", "Astrocyte 2", "Astrocyte 3",
                                     "Dedifferentiating cells", "Fibroblast 1", "Fibroblast 2",
                                     "Senescent cells", "Unknown 1", "Unknown 2"))
seurat@meta.data <- p
Idents(seurat) <- "cell"

## major cell-class grouping, derived from the 18-level subcluster
seurat$cell.major <- factor(seurat$cell)
levels(seurat$cell.major) <- c(rep("Inhibitory neurons", 4), rep("Excitatory neurons", 2),
                                "Dividing cells", "Dividing cells", rep("Glial cells", 4),
                                "Dedifferentiating cells", "Fibroblast", "Fibroblast",
                                "Senescent cells", "Unknown", "Unknown")
seurat$cell.major <- factor(seurat$cell.major,
                             levels = c("Inhibitory neurons", "Excitatory neurons", "Glial cells",
                                        "Dividing cells", "Dedifferentiating cells", "Senescent cells",
                                        "Fibroblast", "Unknown"))

## diagnosis assignment
# unlike the 8-month source, this timepoint's source never explicitly assigns
# a diagnosis column - by the time it's used downstream it's already there,
# so it's probably already in the raw filtered.rds metadata. if it's missing
# here, check the CTRL/sAD donor split in 04_DEG_11month.R and Table 2

save(seurat, file = "/raidixshare_logg01/thais/SarahF/organoids/seurat_11mo_harmony_regressedRibog.rda")

# Figure 2: UMAP + marker heatmaps

df_toPlot <- cbind(Embeddings(object = seurat, reduction = "umap"), seurat@meta.data)
colnames(df_toPlot)[c(1:2)] <- c("UMAP1", "UMAP2")

## color palettes - same major class uses shades of one hue
cols_main <- c(
  "Inhibitory neurons"      = "#3b6b9e",  # blue
  "Excitatory neurons"      = "#6d8c3c",  # olive green
  "Glial cells"             = "#d8b7dd",  # purple
  "Dividing cells"          = "#5abf4b",  # bright green
  "Dedifferentiating cells" = "#e08214",  # orange
  "Senescent cells"         = "#8c510a",  # brown
  "Fibroblast"              = "#b2182b",  # dark red
  "Unknown"                 = "grey80"
)

cols_sub <- c(
  "#1f4e79", "#3b6b9e", "#6a9cc7", "#a1c4e0",  # Inhibitory 1-4 (dark -> pale blue)
  "#3b6b2f", "#6d8c3c", "#8fb35a", "#b5d48a",  # Excitatory 1-4 (dark -> pale green)
  "#7b3294", "#d8b7dd",                        # Glial 1-2 (dark -> pale purple)
  "#5abf4b",                                    # Dividing
  "#e08214", "#fdb863",                        # Dedifferentiating 1-2 (dark -> pale orange)
  "#e06060",                                    # Senescent
  "#b2182b", "#8c510a",                        # Fibroblast 1-2 (dark -> pale red)
  "grey80", "grey70"                            # Unknown 1-2
)

## Fig 2A/2C - UMAP, major class and subcluster
p_sub <- ggplot(df_toPlot, aes(x = UMAP1, y = UMAP2, color = cell)) +
  geom_point(size = 0.6, stroke = 0) +
  scale_color_manual(values = cols_sub) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 10),
    legend.key.size = unit(0.45, "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))

p_main <- ggplot(df_toPlot, aes(x = UMAP1, y = UMAP2, color = cell.major)) +
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

pdf("/raidixshare_logg01/jchen/project/01_sarah_proj/05_organoids/02_results/05_clustering/clustering_umap/11mo_umap.pdf", width = 12, height = 6)
p_sub + p_main
dev.off()

## Fig 2B/2D - marker heatmaps (z-score of mean expression)
marker_genes_major <- c("TUBB3", "MAP2", "GAD1", "GAD2", "SLC17A7", "SLC17A6",
                         "SOX9", "CENPK", "OLIG1", "PDGFRA", "GFAP", "AQP4",
                         "PAX6", "KLF4", "CDKN1A")

marker_genes_sub <- c("TUBB3", "MAP2", "GAD1", "GAD2", "SLC17A7",
                       "SOX9", "CENPK", "OLIG1", "PDGFRA", "GFAP", "AQP4",
                       "PAX6", "CDKN1A", "RPL8")

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

heatmap_major <- make_marker_heatmap(seurat, "cell.major", marker_genes_major)
heatmap_sub   <- make_marker_heatmap(seurat, "cell", marker_genes_sub)

pdf("/raidixshare_logg01/jchen/project/01_sarah_proj/05_organoids/02_results/05_clustering/clustering_umap/11mo_heatmap.pdf", width = 8, height = 8)
grid::grid.newpage()
draw(heatmap_major, padding = unit(c(10, 10, 5, 10), "mm"), newpage = FALSE)
grid::grid.newpage()
draw(heatmap_sub, padding = unit(c(10, 10, 5, 10), "mm"), newpage = FALSE)
dev.off()
