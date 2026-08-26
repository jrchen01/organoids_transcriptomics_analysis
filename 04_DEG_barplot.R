# Number of DEGs per cell subtype
## reads the sig-DEG lists from 04_DEG_8month.R / 04_DEG_11month.R

library(dplyr)
library(ggplot2)
library(tidyr)

## 8-month data (Fig 4A)
mo8_sig_res <- readRDS("/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/04_DEG/8mo_deseq2_sig_DEG.rds")

deg_counts_mo8 <- lapply(names(mo8_sig_res), function(x) {
  celltype  <- sub("_up$|_down$", "", x)
  direction <- ifelse(grepl("_up$", x), "Up", "Down")
  n         <- nrow(mo8_sig_res[[x]])
  data.frame(celltype = celltype, direction = direction, n = n)
})
deg_df_mo8 <- do.call(rbind, deg_counts_mo8)
deg_df_mo8$n <- ifelse(deg_df_mo8$direction == "Down", -deg_df_mo8$n, deg_df_mo8$n)

deg_df_mo8_filtered <- deg_df_mo8 %>% filter(celltype != "Mixed_cells")
deg_df_mo8_filtered$celltype <- factor(deg_df_mo8_filtered$celltype,
  levels = c("Astrocyte", "Dividing_cells", "Excitatory_1",
             "Excitatory_2", "Inhibitory_1", "Inhibitory_2", "OPC"))

cell_colors_mo8 <- c(
  "Astrocyte"      = "#E8A87C",
  "Dividing_cells" = "#C4A35A",
  "Excitatory_1"   = "#6BAF92",
  "Excitatory_2"   = "#4E9B8F",
  "Inhibitory_1"   = "#7BA7D4",
  "Inhibitory_2"   = "#5E8FC4",
  "OPC"            = "#D4A5C9"
)

p_mo8 <- ggplot(deg_df_mo8_filtered, aes(x = n, y = celltype, fill = celltype)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
  geom_text(aes(label = abs(n), hjust = ifelse(n < 0, 1.2, -0.2)), size = 3.5) +
  scale_fill_manual(values = cell_colors_mo8) +
  scale_x_continuous(
    labels = abs,
    name = "Number of DEGs",
    limits = c(-max(abs(deg_df_mo8_filtered$n)) * 1.15, max(abs(deg_df_mo8_filtered$n)) * 1.15),
    expand = expansion(mult = 0)
  ) +
  annotate("text", x = -max(abs(deg_df_mo8_filtered$n)) * 1.1, y = Inf, label = "Down", hjust = -0.2, vjust = 2, size = 3.5, fontface = "bold") +
  annotate("text", x =  max(abs(deg_df_mo8_filtered$n)) * 1.1, y = Inf, label = "Up",   hjust = 1.2,  vjust = 2, size = 3.5, fontface = "bold") +
  labs(y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.y  = element_text(size = 10, face = "bold"),
    axis.text.x  = element_text(size = 9,  face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid   = element_blank()
  )

pdf("/raidixshare_logg01/jchen/project/01_sarah_proj/05_organoids/02_results/04_DEG/figures/number_of_degs_mo8.pdf", width = 8, height = 5)
print(p_mo8)
dev.off()

## 11-month data (Fig 4B)
mo11_sig_res <- readRDS("/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/04_DEG/11mo_deseq2_sig_DEG.rds")

deg_counts_mo11 <- lapply(names(mo11_sig_res), function(x) {
  celltype  <- sub("_up$|_down$", "", x)
  direction <- ifelse(grepl("_up$", x), "Up", "Down")
  n         <- nrow(mo11_sig_res[[x]])
  data.frame(celltype = celltype, direction = direction, n = n)
})
deg_df_mo11 <- do.call(rbind, deg_counts_mo11)
deg_df_mo11$n <- ifelse(deg_df_mo11$direction == "Down", -deg_df_mo11$n, deg_df_mo11$n)

deg_df_mo11_filtered <- deg_df_mo11 %>%
  filter(!celltype %in% c("Dedifferentiating_cells", "Senescent_cells", "Unknown_1"))

deg_df_mo11_filtered$celltype <- factor(deg_df_mo11_filtered$celltype,
  levels = c("Astrocyte_1", "Astrocyte_2", "Astrocyte_3",
             "Dividing_cells", "Dividing_NPC",
             "Excitatory_1", "Excitatory_2",
             "Inhibitory_1", "Inhibitory_2", "Inhibitory_3", "Inhibitory_4",
             "Fibroblast_1", "Fibroblast_2",
             "OPC"))

cell_colors_mo11 <- c(
  "Astrocyte_1"    = "#E8A87C",
  "Astrocyte_2"    = "#D4854A",
  "Astrocyte_3"    = "#C06830",
  "Dividing_cells" = "#C4A35A",
  "Dividing_NPC"   = "#A8873E",
  "Excitatory_1"   = "#6BAF92",
  "Excitatory_2"   = "#4E9B8F",
  "Inhibitory_1"   = "#7BA7D4",
  "Inhibitory_2"   = "#5E8FC4",
  "Inhibitory_3"   = "#4272A8",
  "Inhibitory_4"   = "#2E5A8C",
  "Fibroblast_1"   = "#B07BBF",
  "Fibroblast_2"   = "#8A5A9E",
  "OPC"            = "#D4A5C9"
)

p_mo11 <- ggplot(deg_df_mo11_filtered, aes(x = n, y = celltype, fill = celltype)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
  geom_text(aes(label = abs(n), hjust = ifelse(n < 0, 1.2, -0.2)), size = 3.5) +
  scale_fill_manual(values = cell_colors_mo11) +
  scale_x_continuous(
    labels = abs,
    name = "Number of DEGs",
    limits = c(-max(abs(deg_df_mo11_filtered$n)) * 1.15, max(abs(deg_df_mo11_filtered$n)) * 1.15),
    expand = expansion(mult = 0)
  ) +
  annotate("text", x = -max(abs(deg_df_mo11_filtered$n)) * 1.1, y = Inf, label = "Down", hjust = -0.2, vjust = 2, size = 3.5, fontface = "bold") +
  annotate("text", x =  max(abs(deg_df_mo11_filtered$n)) * 1.1, y = Inf, label = "Up",   hjust = 1.2,  vjust = 2, size = 3.5, fontface = "bold") +
  labs(y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.y  = element_text(size = 10, face = "bold"),
    axis.text.x  = element_text(size = 9,  face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid   = element_blank(),
    plot.margin  = margin(t = 10, r = 10, b = 10, l = 30)
  )

pdf("/raidixshare_logg01/jchen/project/01_sarah_proj/05_organoids/02_results/04_DEG/figures/number_of_degs_mo11.pdf", width = 8, height = 8)
print(p_mo11)
dev.off()
