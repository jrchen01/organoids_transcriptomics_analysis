# SOX9 vs translation module score - 11-month (>10 to 13 months) organoids
# Jiarui Chen

library(Seurat)
library(dplyr)
library(ggplot2)

load("/raidixshare_logg01/thais/SarahF/organoids/seurat_11mo_harmony_regressedRibog.rda")

## Fig 3F - cell-cycle phase UMAP
df_toPlot <- cbind(Embeddings(object = seurat, reduction = "umap"), seurat@meta.data)
colnames(df_toPlot)[c(1:2)] <- c("UMAP1", "UMAP2")

p_phase_mo11 <- ggplot(df_toPlot, aes(x = UMAP1, y = UMAP2, color = Phase)) +
  geom_point(size = 0.6, stroke = 0) +
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

pdf("/raidixshare_logg01/jchen/project/01_sarah_proj/05_organoids/02_results/05_clustering/clustering_umap/11mo_umap_phase.pdf", width = 8, height = 8)
p_phase_mo11
dev.off()

## GOBP translation module score (ribosomal protein genes)
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

## Fig 3E - UMAP of the translation module score
p <- FeaturePlot(seurat, features = "GOBP_Translation1", min.cutoff = "q10", combine = FALSE)

p_gobp_mo11 <- lapply(p, function(x) {
  x +
    labs(title = "GOBP Translation: module score", color = "Module\nScore") +
    theme(
      aspect.ratio = 1,
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line = element_blank(),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text  = element_text(size = 10),
      legend.key.size = unit(0.8, "cm"),
      plot.title = element_text(size = 14, face = "bold")
    )
})

pdf("/raidixshare_logg01/jchen/project/01_sarah_proj/05_organoids/02_results/06_SOX9/11mo_gobp_featureplot.pdf", width = 8, height = 8)
print(p_gobp_mo11)
dev.off()

## Fig 3D - SOX9 expression vs translation module score correlation
a <- FetchData(seurat, vars = c("GOBP_Translation1", "SOX9", "HIF1A"))
cor_result <- cor.test(a$GOBP_Translation1, a$SOX9)
print(cor_result)

p_sox9_cor_mo11 <- ggplot(a, aes(SOX9, GOBP_Translation1)) +
  geom_point() +
  geom_smooth(method = "lm", col = "red") +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 10),
    legend.key.size = unit(0.8, "cm"),
    plot.title = element_text(size = 14, face = "bold")
  ) +
  labs(x = "SOX9 expression", y = "Translation rate",
       title = "Correlation between SOX9 expression and translation") +
  annotate("text", x = 3, y = 4,
           label = sprintf("Pearson correlation coefficient = %.2f\np-value = %.2e",
                            cor_result$estimate, cor_result$p.value))

pdf("/raidixshare_logg01/jchen/project/01_sarah_proj/05_organoids/02_results/06_SOX9/11mo_sox9_corr.pdf", width = 8, height = 8)
p_sox9_cor_mo11
dev.off()
