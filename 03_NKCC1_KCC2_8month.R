# NKCC1 (SLC12A2) / KCC2 (SLC12A5) chloride transporter expression - 8-month organoids


library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)

my18colors <- c(
  "#7B3294", "#2AAD8C", "#3C78D8", "#F4C20D", "#E91E63", "#4CAF50",
  "#FF6F00", "#00BCD4", "#D81B60", "#FF9800", "#5E35B1", "#8D6E63",
  "#AB47BC", "#26A69A", "#1E88E5", "#FDD835", "#EC407A", "#66BB6A"
)

## violin plot of one gene's normalized expression across cell types
plot_gene_violin <- function(seurat_obj, gene, cell_type_col = "cell",
                              gene_label = NULL, colors = my18colors) {
  plot_data <- FetchData(seurat_obj, vars = c(cell_type_col, gene), slot = "data")
  colnames(plot_data) <- c("CellType", "Expression")

  if (is.null(gene_label)) gene_label <- gene

  p <- ggplot(plot_data, aes(x = CellType, y = Expression, fill = CellType)) +
    geom_violin(scale = "width", trim = FALSE, alpha = 0.85, color = "black", linewidth = 0.4) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7, color = "black", linewidth = 0.4, fill = "white") +
    geom_jitter(width = 0.2, size = 0.05, alpha = 0.4, color = "black") +
    scale_fill_manual(values = colors) +
    theme_classic(base_size = 14) +
    guides(fill = guide_legend(override.aes = list(shape = NA, linetype = 0),
                                keywidth = unit(12, "mm"), keyheight = unit(4, "mm"))) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, color = "black", size = 12, face = "bold"),
      axis.text.y = element_text(color = "black", size = 12),
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      legend.position = "none",
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.margin = margin(15, 15, 15, 15)
    ) +
    labs(x = "Cell Type", y = "Normalized Expression", title = gene_label)

  stats <- plot_data %>%
    group_by(CellType) %>%
    summarise(
      mean = mean(Expression),
      Q1 = quantile(Expression, 0.25),
      median = median(Expression),
      Q3 = quantile(Expression, 0.75),
      pct_expressing = sum(Expression > 0) / n() * 100,
      pct_zero = sum(Expression == 0) / n() * 100,
      n_cells = n()
    )

  list(plot = p, data = plot_data, stats = stats)
}

## split by sAD ("AD") vs CTRL ("NC")
plot_gene_by_diagnosis <- function(seurat_obj, gene, cell_type_col = "cell",
                                    diagnosis_col = "diagnosis", gene_label = NULL,
                                    ad_color = "#E64B35", ctrl_color = "#4DBBD5") {
  plot_data <- FetchData(seurat_obj, vars = c(cell_type_col, diagnosis_col, gene), slot = "data")
  colnames(plot_data) <- c("CellType", "Diagnosis", "Expression")

  plot_data$Diagnosis <- dplyr::recode(plot_data$Diagnosis, "CTRL" = "NC", "sAD" = "AD")
  plot_data$Diagnosis <- factor(plot_data$Diagnosis, levels = c("AD", "NC"))

  if (is.null(gene_label)) gene_label <- gene
  diagnosis_colors <- c("AD" = ad_color, "NC" = ctrl_color)

  p <- ggplot(plot_data, aes(x = CellType, y = Expression, fill = Diagnosis)) +
    geom_violin(scale = "width", trim = FALSE, alpha = 0.85, color = "black", linewidth = 0.4,
                position = position_dodge(width = 0.9)) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7, color = "black", linewidth = 0.4,
                 position = position_dodge(width = 0.9)) +
    geom_point(size = 0.05, alpha = 0.4, color = "black",
               position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9)) +
    scale_fill_manual(values = diagnosis_colors) +
    guides(fill = guide_legend(override.aes = list(shape = NA, linetype = 0),
                                keywidth = unit(12, "mm"), keyheight = unit(4, "mm"))) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, color = "black", size = 12, face = "bold"),
      axis.text.y = element_text(color = "black", size = 12),
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.key = element_rect(color = "black", linewidth = 0.3),
      legend.background = element_blank(),
      legend.margin = margin(5, 5, 5, 5),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.margin = margin(15, 15, 15, 15)
    ) +
    labs(x = "Cell Type", y = "Normalized Expression", title = gene_label)

  stats <- plot_data %>%
    group_by(CellType, Diagnosis) %>%
    summarise(
      mean = mean(Expression),
      Q1 = quantile(Expression, 0.25),
      median = median(Expression),
      Q3 = quantile(Expression, 0.75),
      pct_expressing = sum(Expression > 0) / n() * 100,
      pct_zero = sum(Expression == 0) / n() * 100,
      n_cells = n(),
      .groups = "drop"
    )

  list(plot = p, data = plot_data, stats = stats)
}

load("/raidixshare_logg01/thais/SarahF/organoids/seurat_8mo_harmony_regressedRibog.rda")
seurat <- subset(seurat, subset = subcluster != "Unknown")

## all subclusters, no diagnosis split
result_NKCC1 <- plot_gene_violin(seurat, gene = "SLC12A2", gene_label = "NKCC1 (SLC12A2)", cell_type_col = "subcluster")
result_KCC2  <- plot_gene_violin(seurat, gene = "SLC12A5", gene_label = "KCC2 (SLC12A5)",  cell_type_col = "subcluster")

pdf("/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/02_check_genes/plots/8mo_nkcc1_kcc2_subtype.pdf", width = 11, height = 7)
print(result_NKCC1$plot)
print(result_KCC2$plot)
dev.off()

write.csv(result_NKCC1$stats, row.names = FALSE, file = "/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/02_check_genes/plots/8mo_NKCC1_expression_stats.csv")
write.csv(result_KCC2$stats,  row.names = FALSE, file = "/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/02_check_genes/plots/8mo_KCC2_expression_stats.csv")

## Excitatory/Inhibitory neuron subclusters, sAD vs CTRL
seurat_neuron <- subset(seurat, subset = subcluster %in%
                           c("Excitatory 1", "Excitatory 2", "Inhibitory 1", "Inhibitory 2"))

result_nkcc1_diagnosis <- plot_gene_by_diagnosis(seurat_neuron, gene = "SLC12A2", gene_label = "NKCC1 (SLC12A2)",
                                                  cell_type_col = "subcluster", diagnosis_col = "diagnosis")
result_kcc2_diagnosis  <- plot_gene_by_diagnosis(seurat_neuron, gene = "SLC12A5", gene_label = "KCC2 (SLC12A5)",
                                                  cell_type_col = "subcluster", diagnosis_col = "diagnosis")

pdf("/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/02_check_genes/plots/8mo_nkcc1_kcc2_subtype_diagnosis.pdf", width = 10, height = 7)
print(result_nkcc1_diagnosis$plot)
print(result_kcc2_diagnosis$plot)
dev.off()

write.csv(result_nkcc1_diagnosis$stats, row.names = FALSE, file = "/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/02_check_genes/plots/8mo_NKCC1_expression_diagnosis_stats.csv")
write.csv(result_kcc2_diagnosis$stats,  row.names = FALSE, file = "/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/02_check_genes/plots/8mo_KCC2_expression_diagnosis_stats.csv")

## mean SLC12A2:SLC12A5 ratio by diagnosis, Excitatory 1/2
mean_expr_summary <- result_nkcc1_diagnosis$stats %>%
  rename(mean_SLC12A2 = mean) %>%
  select(CellType, Diagnosis, mean_SLC12A2) %>%
  inner_join(
    result_kcc2_diagnosis$stats %>% rename(mean_SLC12A5 = mean) %>% select(CellType, Diagnosis, mean_SLC12A5),
    by = c("CellType", "Diagnosis")
  ) %>%
  filter(CellType %in% c("Excitatory 1", "Excitatory 2")) %>%
  mutate(
    Diagnosis = recode(Diagnosis, "AD" = "sAD", "NC" = "CTRL"),
    `SLC12A2:SLC12A5` = mean_SLC12A2 / mean_SLC12A5
  )

print(mean_expr_summary)
## this whole block isn't from the original notebooks - wrote it to reproduce
## the "Mean Expression" table in the report; path follows the other 8mo outputs
write.csv(mean_expr_summary, row.names = FALSE, file = "/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/02_check_genes/plots/8mo_NKCC1_KCC2_mean_expression_ratio_summary.csv")

rm(seurat, seurat_neuron)

# per-patient pseudobulk SLC12A2:SLC12A5 ratio - Excitatory 1/2 only
load("/raidixshare_logg01/thais/SarahF/organoids/seurat_8mo_harmony_regressedRibog.rda")

exc <- subset(seurat, subset = subcluster %in% c("Excitatory 1", "Excitatory 2"))
DefaultAssay(exc) <- "RNA"

expr <- GetAssayData(exc, slot = "data")[c("SLC12A2", "SLC12A5"), ]
expr <- as.data.frame(t(as.matrix(expr)))
df <- cbind(exc@meta.data, expr)

pseudo <- df %>%
  group_by(patient_id, diagnosis, subcluster) %>%
  summarise(
    n_cells = n(),
    mean_SLC12A2 = mean(SLC12A2),
    mean_SLC12A5 = mean(SLC12A5),
    .groups = "drop"
  ) %>%
  filter(n_cells >= 10) %>%
  mutate(
    log_ratio = log2(mean_SLC12A2 + 1e-3) - log2(mean_SLC12A5 + 1e-3),
    ratio = (mean_SLC12A2 + 1e-3) / (mean_SLC12A5 + 1e-3)
  )

boxplot_theme <- theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4)
  )

ex1 <- subset(pseudo, subcluster == "Excitatory 1")
# print(wilcox.test(log_ratio ~ diagnosis, data = ex1))

p_exc1 <- ggplot(ex1, aes(x = diagnosis, y = log_ratio, fill = diagnosis)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.6, linewidth = 0.6) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("grey40", "firebrick")) +
  labs(title = "Excitatory 1", y = "log(SLC12A2/SLC12A5)", x = NULL) +
  boxplot_theme

ex2 <- subset(pseudo, subcluster == "Excitatory 2")
# print(wilcox.test(log_ratio ~ diagnosis, data = ex2))

p_exc2 <- ggplot(ex2, aes(x = diagnosis, y = log_ratio, fill = diagnosis)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.6, linewidth = 0.6) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("grey40", "firebrick")) +
  labs(title = "Excitatory 2", y = "log(SLC12A2/SLC12A5)", x = NULL) +
  boxplot_theme

pdf("/raidixshare_logg01/jchen/project/01_sarah_proj/05_organoids/02_results/02_check_genes/comparing_NKCC1_KCC2_boxplots/8mo_NKCC1_KCC2_EXC_boxplot.pdf", width = 10, height = 4)
p_exc1 + p_exc2
dev.off()
