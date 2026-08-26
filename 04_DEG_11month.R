# Pseudobulk DESeq2 (sAD vs CTRL) - 11-month (>10 to 13 months) organoids

library(Seurat)
library(dplyr)
library(SingleCellExperiment)
library(DESeq2)
library(tidyverse)
library(Matrix.utils)
library(stringr)
library(kableExtra)

source("functions_deseq2_pseudobulk.R")

csv_dir <- "/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/04_DEG/11mo_csv/"
genes_dir <- "/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/04_DEG/11mo_genes/"
sig_deg_rds <- "/home/jchen/data1/jchen/project/01_sarah_proj/05_organoids/02_results/04_DEG/11mo_deseq2_sig_DEG.rds"

load("/raidixshare_logg01/thais/SarahF/organoids/seurat_11mo_harmony_regressedRibog.rda")

## donor sex / APOE genotype, not in the source object's metadata
## only 8 of the 11 donors from 8 months show up here - 3483, ADRC_2800, and the MCI donor 3341 are gone
seurat@meta.data$sex <- NA
seurat@meta.data$apoe <- NA

donor_info <- list(
  "UCI40" = c(sex = "M", apoe = "23"),
  "40"    = c(sex = "M", apoe = "33"),
  "3551"  = c(sex = "F", apoe = "34"),
  "UCI52" = c(sex = "F", apoe = "23"),
  "UCI2"  = c(sex = "M", apoe = "33"),
  "UCI22" = c(sex = "F", apoe = "33"),
  "2991"  = c(sex = "M", apoe = "33"),
  "UCI96" = c(sex = "F", apoe = "44")
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

## collapse "cell" (subcluster) into DESeq2-safe labels, spaces -> underscores
## labels not listed here pass through unchanged (e.g. "OPC")
seurat@meta.data <- seurat@meta.data %>%
  mutate(
    subcluster2 = case_when(
      cell %in% c("Fibroblast 1")           ~ "Fibroblast_1",
      cell %in% c("Astrocyte 1")            ~ "Astrocyte_1",
      cell %in% c("Astrocyte 3")            ~ "Astrocyte_3",
      cell %in% c("Unknown 2")              ~ "Unknown_2",
      cell %in% c("Inhibitory 4")           ~ "Inhibitory_4",
      cell %in% c("Excitatory 1")           ~ "Excitatory_1",
      cell %in% c("Dedifferentiating cells")~ "Dedifferentiating_cells",
      cell %in% c("Dividing NPC")           ~ "Dividing_NPC",
      cell %in% c("Inhibitory 1")           ~ "Inhibitory_1",
      cell %in% c("Astrocyte 2")            ~ "Astrocyte_2",
      cell %in% c("Senescent cells")        ~ "Senescent_cells",
      cell %in% c("Fibroblast 2")           ~ "Fibroblast_2",
      cell %in% c("Dividing cells")         ~ "Dividing_cells",
      cell %in% c("Unknown 1")              ~ "Unknown_1",
      cell %in% c("Inhibitory 3")           ~ "Inhibitory_3",
      cell %in% c("Inhibitory 2")           ~ "Inhibitory_2",
      cell %in% c("Excitatory 2")           ~ "Excitatory_2",
      TRUE ~ cell
    )
  )

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
## Unknown_2 fails here - model matrix not full rank, 1 donor per group

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
  write.table(genes, file = file.path(genes_dir, paste0("11mo_", nm, ".txt")),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}
