# data assessment

library(Seurat)
library(ggplot2)

astro.rna <- readRDS("/raidixshare_logg01/thais/SarahF/organoids/seurat_noMatrigel_0.05.rds")
head(astro.rna@meta.data)

## initial checks
head(astro.rna@meta.data$patient_id)
head(astro.rna@meta.data$diagnosis)
head(astro.rna@meta.data$Time)
head(astro.rna@meta.data$Phase)
head(astro.rna@meta.data$apoe)

## table overview
table(astro.rna@meta.data$patient_id)
astro_rna_cell_counts <- as.data.frame(table(astro.rna@meta.data$patient_id, astro.rna@meta.data$Time))
colnames(astro_rna_cell_counts) <- c("patient_id", "Time", "Cell_Count")
patient_info <- unique(astro.rna@meta.data[, c("patient_id", "diagnosis", "apoe")])
astro_rna_cell_counts <- merge(astro_rna_cell_counts, patient_info, by = "patient_id", all.x = TRUE)
print(astro_rna_cell_counts)

## barplot of counts for each cell line (patient)
patient_data <- as.data.frame(table(astro.rna@meta.data$patient_id))
ggplot(patient_data, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Bar Plot of Patient Counts", x = "Patient ID", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

table(astro.rna@meta.data$patient_id, astro.rna@meta.data$Time)

## barplot of counts for each patient by different time points
patient_time_data <- as.data.frame(table(astro.rna@meta.data$patient_id, astro.rna@meta.data$Time))
patient_time_data <- patient_time_data[patient_time_data$Freq > 0, ]
ggplot(patient_time_data, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = Freq), position = position_dodge(width = 0.9), vjust = -0.5) +
  theme_minimal() +
  labs(title = "Patient ID Counts by Time", x = "Patient ID", y = "Count", fill = "Time") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 