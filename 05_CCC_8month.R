# Cell-cell communication (CellChat) - 8-month organoids
## CellChat run separately per diagnosis group, restricted to "Non-protein Signaling" in CellChatDB.human

library(Seurat)
library(CellChat)
library(dplyr)
library(tidyr)
library(ggalluvial)
library(presto)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(glue)

source("functions_cellchat_plots.R")

# build CellChat objects
load("/raidixshare_logg01/thais/SarahF/organoids/seurat_8mo_harmony_regressedRibog.rda")
seurat <- subset(seurat, subset = !(subcluster %in% c("Unknown", "Mixed cells")))

run_cellchat <- function(seurat_group) {
  matrix_data <- seurat_group@assays$RNA$data %>% na.omit()
  cellchat <- createCellChat(object = matrix_data, meta = seurat_group@meta.data, group.by = "subcluster")
  CellChatDB.use <- subsetDB(CellChatDB.human, search = "Non-protein Signaling", key = "annotation")
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  future::plan("multisession", workers = 4)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat
}

## subtype level (Fig 5a-b)
cellchat_sub_ad <- run_cellchat(subset(seurat, subset = diagnosis == "sAD"))
cellchat_sub_nc <- run_cellchat(subset(seurat, subset = diagnosis == "CTRL"))

saveRDS(cellchat_sub_ad, "/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/01_ccc/ccc_8mo_non_protein_AD_subtype.rds")
saveRDS(cellchat_sub_nc, "/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/01_ccc/ccc_8mo_non_protein_NC_subtype.rds")

## major cell-class level - not in Fig 5, kept for reference
run_cellchat_major <- function(seurat_group) {
  matrix_data <- seurat_group@assays$RNA$data %>% na.omit()
  cellchat <- createCellChat(object = matrix_data, meta = seurat_group@meta.data, group.by = "cell")
  CellChatDB.use <- subsetDB(CellChatDB.human, search = "Non-protein Signaling", key = "annotation")
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  future::plan("multisession", workers = 4)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat
}

cellchat_major_ad <- run_cellchat_major(subset(seurat, subset = diagnosis == "sAD"))
cellchat_major_nc <- run_cellchat_major(subset(seurat, subset = diagnosis == "CTRL"))

saveRDS(cellchat_major_ad, "/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/01_ccc/ccc_8mo_non_protein_AD_major.rds")
saveRDS(cellchat_major_nc, "/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/01_ccc/ccc_8mo_non_protein_NC_major.rds")

# Fig 5a - circle plots, Fig 5b - outgoing/incoming signal bar graphs
ad_data <- extract_signaling_data(cellchat_sub_ad, "AD")
nc_data <- extract_signaling_data(cellchat_sub_nc, "NC")
combined_data <- rbind(ad_data, nc_data)

pdf("/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/01_ccc/plots/8mo_ccc_non_protein_signals_subtype.pdf", width = 12, height = 10)

par(mfrow = c(2, 2), xpd = TRUE)
groupSize_ad <- as.numeric(table(cellchat_sub_ad@idents))
invisible(netVisual_circle(cellchat_sub_ad@net$count, vertex.weight = groupSize_ad, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions (AD)"))
invisible(netVisual_circle(cellchat_sub_ad@net$weight, vertex.weight = groupSize_ad, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength (AD)"))

groupSize_nc <- as.numeric(table(cellchat_sub_nc@idents))
invisible(netVisual_circle(cellchat_sub_nc@net$count, vertex.weight = groupSize_nc, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions (NC)"))
invisible(netVisual_circle(cellchat_sub_nc@net$weight, vertex.weight = groupSize_nc, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength (NC)"))

plot_outgoing_count    <- create_cellchat_academic_plot(combined_data, "count", "outgoing")
plot_incoming_count    <- create_cellchat_academic_plot(combined_data, "count", "incoming")
plot_outgoing_strength <- create_cellchat_academic_plot(combined_data, "strength", "outgoing")
plot_incoming_strength <- create_cellchat_academic_plot(combined_data, "strength", "incoming")

grid.arrange(plot_outgoing_count, plot_incoming_count,
             plot_outgoing_strength, plot_incoming_strength,
             ncol = 2, nrow = 2)
dev.off()

## final combined differential (log2FC) + raw signal-strength bar plot
control_send_weight    <- rowSums(cellchat_sub_nc@net$weight)
control_receive_weight <- colSums(cellchat_sub_nc@net$weight)
ad_send_weight    <- rowSums(cellchat_sub_ad@net$weight)
ad_receive_weight <- colSums(cellchat_sub_ad@net$weight)

all_cell_types <- unique(c(names(control_send_weight), names(ad_send_weight)))

df_send_weight <- data.frame(
  CellType = all_cell_types,
  Control = sapply(all_cell_types, function(x) ifelse(x %in% names(control_send_weight), control_send_weight[x], 0)),
  AD = sapply(all_cell_types, function(x) ifelse(x %in% names(ad_send_weight), ad_send_weight[x], 0))
)
df_receive_weight <- data.frame(
  CellType = all_cell_types,
  Control = sapply(all_cell_types, function(x) ifelse(x %in% names(control_receive_weight), control_receive_weight[x], 0)),
  AD = sapply(all_cell_types, function(x) ifelse(x %in% names(ad_receive_weight), ad_receive_weight[x], 0))
)

ad_diff_send    <- plot_diff_barplot_2group(df_send_weight, all_cell_types, condition = "AD", ref = "Control", threshold = 0.5)
ad_diff_receive <- plot_diff_barplot_2group(df_receive_weight, all_cell_types, condition = "AD", ref = "Control", threshold = 0.5, title = "Differential input (log2FC)")
ad_signal_send    <- plot_signal_barplot_2group(df_send_weight, all_cell_types, condition = "AD", ref = "Control", xlab = "Total Outgoing Signal Strength")
ad_signal_receive <- plot_signal_barplot_2group(df_receive_weight, all_cell_types, condition = "AD", ref = "Control", xlab = "Total Incoming Signal Strength")

combined_ad_send    <- ad_diff_send + ad_signal_send
combined_ad_receive <- ad_diff_receive + ad_signal_receive

pdf("/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/01_ccc/plots/8mo/8mo_ccc_non_protein_signals_subtype_barplot_final.pdf", width = 10, height = 5)
print(combined_ad_send)
print(combined_ad_receive)
dev.off()

# differential chord diagrams (AD - NC), normalized to equal total signal
weight_ad <- cellchat_sub_ad@net$weight
weight_nc <- cellchat_sub_nc@net$weight

if (!all(rownames(weight_ad) == rownames(weight_nc))) {
  stop("Cell types don't match between AD and NC!")
}

norm_factor <- sum(weight_nc) / sum(weight_ad)
diff_matrix <- weight_ad * norm_factor - weight_nc

pdf("/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/01_ccc/plots/8mo/8mo_EX1_outgoing_differential.pdf", width = 8, height = 8)
my_netVisual_chord(
  diff_matrix, sources.use = "Excitatory 1", targets.use = NULL,
  edge.transparency = TRUE, show.legend = TRUE, legend.pos.y = 20,
  title.name = "Excitatory 1 Differential Outgoing Signals\n(AD - NC)",
  vertex.label.cex = 0.8, small.gap = 1, big.gap = 10,
  remove.isolate = FALSE, directional = 1, title.cex = 1.2, range_value = 0.5
)
dev.off()

pdf("/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/01_ccc/plots/8mo/8mo_INH2_outgoing_differential.pdf", width = 8, height = 8)
my_netVisual_chord(
  diff_matrix, sources.use = "Inhibitory 2", targets.use = NULL,
  edge.transparency = TRUE, show.legend = TRUE, legend.pos.y = 20,
  title.name = "Inhibitory 2 Differential Outgoing Signals\n(AD - NC)",
  vertex.label.cex = 0.8, small.gap = 1, big.gap = 10,
  remove.isolate = FALSE, directional = 1, title.cex = 1.2, range_value = 0.5
)
dev.off()
