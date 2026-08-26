# CellChat supplementary - LR pair / pathway level differential signaling, 11-month organoids
## GABA-B and Glutamate pathways for Inhibitory 3 and Inhibitory 4
## input: subtype-level CellChat objects from 05_CCC_11month.R

library(Seurat)
library(CellChat)
library(dplyr)
library(ggplot2)
library(glue)

source("functions_cellchat_plots.R")

cellchat_sub_ad <- readRDS("/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/01_ccc/ccc_11mo_non_protein_AD_subtype.rds")
cellchat_sub_nc <- readRDS("/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/01_ccc/ccc_11mo_non_protein_NC_subtype.rds")

## condition 1 = reference (NC), condition 2 = test (AD)
cellchat_1 <- cellchat_sub_nc
cellchat_2 <- cellchat_sub_ad

## aggregate one LR pair or pathway into a normalized AD-minus-NC difference matrix
aggregate_pairs <- function(pairLR.use.name, agg_method = "weight_threshold", cut_off = 0.01) {
  net1 <- cellchat_1@net$prob
  net2 <- cellchat_2@net$prob
  global_norm_factor <- sum(cellchat_1@net$prob) / sum(cellchat_2@net$prob)

  net1_lr <- net1[, , pairLR.use.name]
  net2_lr <- net2[, , pairLR.use.name]

  if (agg_method == "weight_threshold") {
    net1_lr[net1_lr < cut_off] <- 0
    net2_lr[net2_lr < cut_off] <- 0
  }

  net2_lr * global_norm_factor - net1_lr
}

aggregate_pairs_pathway <- function(pairLR.use.name, agg_method = "weight_threshold", cut_off = 0.01) {
  net1 <- cellchat_1@netP$prob
  net2 <- cellchat_2@netP$prob
  global_norm_factor <- sum(cellchat_1@net$prob) / sum(cellchat_2@net$prob)

  net1_lr <- net1[, , pairLR.use.name]
  net2_lr <- net2[, , pairLR.use.name]

  if (agg_method == "weight_threshold") {
    net1_lr[net1_lr < cut_off] <- 0
    net2_lr[net2_lr < cut_off] <- 0
  }

  net2_lr * global_norm_factor - net1_lr
}

plot_and_save_pathway <- function(source_celltype, pathway, save_path) {
  net.diff <- aggregate_pairs_pathway(pairLR.use.name = pathway)
  plt_chord <- my_netVisual_chord(net.diff, sources.use = source_celltype,
                                   edge.transparency = TRUE, show.legend = TRUE,
                                   title.name = sprintf("%s\n%s Pathway (AD - NC)", source_celltype, pathway))
  pdf(glue("{save_path}/{pathway}.pdf"), width = 8, height = 8)
  print(plt_chord)
  dev.off()
}

plot_and_save_lr_pair <- function(source_celltype, lr_pair, save_path) {
  net.diff <- aggregate_pairs(pairLR.use.name = lr_pair)
  plt_chord <- my_netVisual_chord(net.diff, sources.use = source_celltype,
                                   edge.transparency = TRUE, show.legend = TRUE,
                                   title.name = lr_pair)
  pdf(glue("{save_path}/{lr_pair}.pdf"), width = 8, height = 8)
  print(plt_chord)
  dev.off()
}

## Inhibitory 4 - GABA-B and Glutamate pathway-level differential signaling
inh4_path <- "/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/01_ccc/plots/11mo/INH4_LR"
plot_and_save_pathway("Inhibitory 4", "GABA-B", inh4_path)
plot_and_save_pathway("Inhibitory 4", "Glutamate", inh4_path)

## Inhibitory 4 - Glutamate LR pairs ranked by AD-NC difference
net_ad <- cellchat_sub_ad@net$prob
net_nc <- cellchat_sub_nc@net$prob
all_LR <- dimnames(net_ad)[[3]]
glu_LR <- all_LR[grepl("Glutamate|GLU|GRIA|GRIN|GRM", all_LR, ignore.case = TRUE)]

df_glu_inh4 <- data.frame()
for (lr in glu_LR) {
  ad_val <- sum(net_ad["Inhibitory 4", , lr], na.rm = TRUE)
  nc_val <- sum(net_nc["Inhibitory 4", , lr], na.rm = TRUE)
  df_glu_inh4 <- rbind(df_glu_inh4, data.frame(LR_pair = lr, AD = ad_val, NC = nc_val,
                                                Diff = ad_val - nc_val, AbsDiff = abs(ad_val - nc_val)))
}
df_glu_inh4 <- df_glu_inh4[order(-df_glu_inh4$AbsDiff), ]
print(df_glu_inh4)

glu_lr_pairs_inh4 <- c("Glutamate-Glu-SLC1A1_GLS_GRIA2",
                        "Glutamate-Glu-SLC1A3_GLS_GRIA2",
                        "Glutamate-Glu-SLC1A2_GLS_GRIA2")
for (lr in glu_lr_pairs_inh4) plot_and_save_lr_pair("Inhibitory 4", lr, inh4_path)

## Inhibitory 4 - GABA LR pairs ranked by AD-NC difference
gaba_LR <- all_LR[grepl("GABA", all_LR, ignore.case = TRUE)]

df_gaba_inh4 <- data.frame()
for (lr in gaba_LR) {
  ad_val <- sum(net_ad["Inhibitory 4", , lr], na.rm = TRUE)
  nc_val <- sum(net_nc["Inhibitory 4", , lr], na.rm = TRUE)
  df_gaba_inh4 <- rbind(df_gaba_inh4, data.frame(LR_pair = lr, AD = ad_val, NC = nc_val,
                                                  Diff = ad_val - nc_val, AbsDiff = abs(ad_val - nc_val)))
}
df_gaba_inh4 <- df_gaba_inh4[order(-df_gaba_inh4$AbsDiff), ]
print(df_gaba_inh4)

gaba_lr_pairs_inh4 <- c("GABA-B-GAD1_SLC6A1_GABBR1", "GAD1_SLC6A1_GABBR2")
for (lr in gaba_lr_pairs_inh4) plot_and_save_lr_pair("Inhibitory 4", lr, inh4_path)

## Inhibitory 3 - Glutamate pathway-level and LR pair signaling
inh3_path <- "/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/01_ccc/plots/11mo/INH3_LR"
plot_and_save_pathway("Inhibitory 3", "Glutamate", inh3_path)

df_glu_inh3 <- data.frame()
for (lr in glu_LR) {
  ad_val <- sum(net_ad["Inhibitory 3", , lr], na.rm = TRUE)
  nc_val <- sum(net_nc["Inhibitory 3", , lr], na.rm = TRUE)
  df_glu_inh3 <- rbind(df_glu_inh3, data.frame(LR_pair = lr, AD = ad_val, NC = nc_val,
                                                Diff = ad_val - nc_val, AbsDiff = abs(ad_val - nc_val)))
}
df_glu_inh3 <- df_glu_inh3[order(-df_glu_inh3$AbsDiff), ]
print(df_glu_inh3)

glu_lr_pairs_inh3 <- c("Glutamate-Glu-SLC1A3_GLS_GRIA2", "Glutamate-Glu-SLC1A3_GLS_GRIK2")
for (lr in glu_lr_pairs_inh3) plot_and_save_lr_pair("Inhibitory 3", lr, inh3_path)
