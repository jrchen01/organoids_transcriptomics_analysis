# Pseudobulk DESeq2 (sAD vs CTRL) - 8-month organoids

library(Seurat)
library(dplyr)
library(SingleCellExperiment)
library(DESeq2)
library(tidyverse)
library(Matrix.utils)
library(stringr)
library(kableExtra)

source("functions_deseq2_pseudobulk.R")

csv_dir <- "/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/04_DEG/8mo_csv/"
genes_dir <- "/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/04_DEG/8mo_genes/"
sig_deg_rds <- "/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/04_DEG/8mo_deseq2_sig_DEG.rds"

load("/raidixshare_logg01/thais/SarahF/organoids/seurat_8mo_harmony_regressedRibog.rda")

## collapse subcluster into DEG labels (merge OPC 1/2 -> OPC, spaces -> underscores)
seurat@meta.data <- seurat@meta.data %>%
  mutate(
    subcluster2 = case_when(
      subcluster %in% c("OPC 1", "OPC 2")   ~ "OPC",
      subcluster %in% c("Dividing cells")   ~ "Dividing_cells",
      subcluster %in% c("Inhibitory 1")     ~ "Inhibitory_1",
      subcluster %in% c("Inhibitory 2")     ~ "Inhibitory_2",
      subcluster %in% c("Excitatory 1")     ~ "Excitatory_1",
      subcluster %in% c("Excitatory 2")     ~ "Excitatory_2",
      subcluster %in% c("Mixed cells")      ~ "Mixed_cells",
      TRUE ~ subcluster
    )
  )

## donor sex / APOE genotype, not in the source object's metadata
seurat@meta.data$sex <- NA
seurat@meta.data$apoe <- NA

donor_info <- list(
  "UCI40"     = c(sex = "M", apoe = "23"),
  "40"        = c(sex = "M", apoe = "33"),
  "3551"      = c(sex = "F", apoe = "34"),
  "3483"      = c(sex = "F", apoe = "34"),
  "UCI52"     = c(sex = "F", apoe = "23"),
  "UCI2"      = c(sex = "M", apoe = "33"),
  "UCI22"     = c(sex = "F", apoe = "33"),
  "ADRC_2800" = c(sex = "M", apoe = "34"),
  "2991"      = c(sex = "M", apoe = "33"),
  "UCI96"     = c(sex = "F", apoe = "44"),
  "3341"      = c(sex = "M", apoe = "34")  # MCI donor, excluded downstream (diagnosis %in% c("CTRL","sAD") only)
)
for (pid in names(donor_info)) {
  seurat@meta.data[seurat@meta.data$patient_id == pid, "sex"]  <- donor_info[[pid]]["sex"]
  seurat@meta.data[seurat@meta.data$patient_id == pid, "apoe"] <- donor_info[[pid]]["apoe"]
}

table_data <- seurat@meta.data %>%
  group_by(patient_id, diagnosis, sex, apoe) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  arrange(diagnosis, patient_id)
print(table_data)

## run pseudobulk DESeq2 per cell subtype
all_celltypes <- unique(seurat$subcluster2)
celltypes_to_analyze <- all_celltypes[all_celltypes != "Unknown"]
cat("Cell types to analyze:\n")
print(celltypes_to_analyze)

de_results_list <- lapply(celltypes_to_analyze, function(ct) {
  cat("Processing cell type:", ct, "\n")
  tryCatch({
    result <- run_pseudobulk_deseq2(
      seurat_obj = seurat,
      celltype_label = "subcluster2",
      celltype = ct,
      diagnosis_label = "diagnosis",
      diagnosis_group = "sAD",
      reference_group = "CTRL",
      patient_id_label = "patient_id",
      OUTDIR = csv_dir
    )
    cat("Significant genes found:", result$n_significant, "\n\n")
    result
  }, error = function(e) {
    cat("ERROR in cell type", ct, ":", e$message, "\n\n")
    NULL
  })
})
names(de_results_list) <- paste0("de_result_", celltypes_to_analyze)
## Fibroblast fails here - model matrix not full rank, only 5 donors (2 vs 3)

## extract significant genes: Wald test, pvalue < 0.05, split by direction
significant_results <- list()

for (result_name in names(de_results_list)) {
  celltype <- de_results_list[[result_name]]$celltype
  all_res <- de_results_list[[result_name]]$all_results
  if (is.null(all_res)) {
    message("Skipping ", result_name, ": all_results is NULL")
    next
  }

  wald_results <- all_res %>% filter(DESeq_test == "Wald" & pvalue < 0.05)

  significant_results[[paste0(celltype, "_up")]]   <- wald_results %>% filter(log2FoldChange > 0.25)
  significant_results[[paste0(celltype, "_down")]] <- wald_results %>% filter(log2FoldChange < -0.25)
}

saveRDS(significant_results, sig_deg_rds)

## gene lists per celltype/direction, for GO enrichment in Metascape
for (nm in names(significant_results)) {
  genes <- unique(significant_results[[nm]]$gene)
  write.table(genes, file = file.path(genes_dir, paste0("8mo_", nm, ".txt")),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}
