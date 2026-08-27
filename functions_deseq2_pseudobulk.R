# pseudobulk DESeq2 helper


library(Seurat)
library(SingleCellExperiment)
library(DESeq2)
library(tidyverse)
library(Matrix.utils)
library(stringr)

run_pseudobulk_deseq2 <- function(seurat_obj,
                                   celltype_label = "celltype_label",
                                   celltype = "EX_1",
                                   diagnosis_label = "diagnosis",
                                   diagnosis_group = "AD",
                                   reference_group = "NC",
                                   patient_id_label = "patient_id",
                                   log2FC = 0.25,
                                   padj_cutoff = 0.05,
                                   min_cells_per_donor = 10,
                                   min_count_threshold = 10,
                                   min_samples_with_count = 2,
                                   run_both_tests = TRUE,
                                   OUTDIR = "/path/to/output/") {

  prepare_dds <- function(seurat, cell_type) {
    cat("Processing cell type:", cell_type, "\n")

    p <- seurat@meta.data
    p <- p[p[[celltype_label]] %in% c(cell_type), ]

    if (nrow(p) == 0) {
      stop("No cells found for cell type: ", cell_type)
    }

    required_cols <- c(patient_id_label, diagnosis_label)
    missing_cols <- setdiff(required_cols, colnames(p))
    if (length(missing_cols) > 0) {
      stop("Missing required columns in metadata: ", paste(missing_cols, collapse = ", "))
    }

    # keep only donors contributing enough cells to this cell type
    ids <- table(p[[patient_id_label]])
    ids <- names(ids[ids >= min_cells_per_donor])

    if (length(ids) < 2) {
      stop("Not enough donors with sufficient cells (>= ", min_cells_per_donor, " cells)")
    }

    p <- p[p[[patient_id_label]] %in% ids, ]
    cat("Selected", length(ids), "donors with >=", min_cells_per_donor, "cells\n")

    # raw counts for this cell type's cells
    counts <- seurat@assays$RNA$counts
    counts <- counts[, rownames(p)]

    sce <- SingleCellExperiment(assays = list(counts = counts),
                                 colData = p[, patient_id_label, drop = FALSE])

    # pseudobulk: sum raw counts across all cells of this type, per donor
    groups <- SingleCellExperiment::colData(sce)
    pb <- Matrix.utils::aggregate.Matrix(t(counts(sce)),
                                          groupings = groups, fun = "sum")

    df <- as.data.frame(t(pb))

    # one metadata row per donor
    p <- p[!duplicated(p[[patient_id_label]]), ]
    rownames(p) <- p[[patient_id_label]]

    p[[diagnosis_label]] <- factor(p[[diagnosis_label]])
    p[[diagnosis_label]] <- relevel(p[[diagnosis_label]], ref = reference_group)

    p <- p[colnames(df), ]

    if (!identical(rownames(p), colnames(df))) {
      stop("Error: rownames of p and colnames of df do not match!")
    }
    if (nrow(p) != ncol(df)) {
      stop("Error: Number of rows in p does not match number of columns in df!")
    }

    cat("Sample distribution:\n")
    print(table(p[[diagnosis_label]]))

    return(list(df = df, p = p))
  }

  result_data <- prepare_dds(seurat_obj, celltype)

  # design: ~ sex + diagnosis, to control for gender-related confounding
  dds <- DESeqDataSetFromMatrix(result_data$df,
                                 colData = result_data$p,
                                 design = as.formula(paste("~sex + ", diagnosis_label)))

  # keep genes with >= min_count_threshold counts in >= min_samples_with_count samples
  keep <- rowSums(counts(dds) >= min_count_threshold) >= min_samples_with_count
  dds <- dds[keep, ]

  message(str_glue("-------------------------------------------------------\n\nRunning DESeq2: {diagnosis_group} vs {reference_group} in {celltype}
                   \n------------------------------------------------------"))

  results_list <- list()

  if (run_both_tests) {
    cat("Running LRT test...\n")
    dds_LRT <- DESeq(dds, test = "LRT", reduced = as.formula("~sex"))
    res_LRT <- results(dds_LRT, contrast = c(diagnosis_label, diagnosis_group, reference_group))

    res_df_LRT <- as.data.frame(res_LRT) %>%
      rownames_to_column(var = "gene") %>%
      as_tibble() %>%
      arrange(padj, pvalue) %>%
      mutate(
        DESeq_test = "LRT",
        celltype = celltype,
        diagnosis = diagnosis_group,
        comparison = paste(diagnosis_group, "vs", reference_group)
      )

    results_list[["LRT"]] <- res_df_LRT
  }

  cat("Running Wald test...\n")
  dds_Wald <- DESeq(dds, test = "Wald")
  res_Wald <- results(dds_Wald, contrast = c(diagnosis_label, diagnosis_group, reference_group))

  res_df_Wald <- as.data.frame(res_Wald) %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    arrange(padj, pvalue) %>%
    mutate(
      DESeq_test = "Wald",
      celltype = celltype,
      diagnosis = diagnosis_group,
      comparison = paste(diagnosis_group, "vs", reference_group)
    )

  results_list[["Wald"]] <- res_df_Wald

  if (run_both_tests) {
    res_df_all <- bind_rows(results_list$LRT, results_list$Wald)
  } else {
    res_df_all <- results_list$Wald
  }

  res_df_all <- res_df_all %>%
    mutate(
      significant = padj < padj_cutoff & abs(log2FoldChange) >= log2FC,
      direction = case_when(
        padj < padj_cutoff & log2FoldChange >= log2FC ~ "Up",
        padj < padj_cutoff & log2FoldChange <= -log2FC ~ "Down",
        TRUE ~ "NS"
      )
    )

  res_df_all_sig <- res_df_all %>%
    filter(padj < padj_cutoff & abs(log2FoldChange) > log2FC)

  timestamp <- format(Sys.time(), "%y%m%d")
  all_results_file <- file.path(OUTDIR, paste0(celltype, "_", diagnosis_group, "_vs_", reference_group, "_all_results_", timestamp, ".csv"))
  sig_results_file <- file.path(OUTDIR, paste0(celltype, "_", diagnosis_group, "_vs_", reference_group, "_sig_results_", timestamp, ".csv"))

  write_csv(res_df_all, all_results_file)
  write_csv(res_df_all_sig, sig_results_file)

  cat("\n=== Results Summary ===\n")
  cat("Total genes tested:", length(unique(res_df_all$gene)), "\n")
  cat("Significant genes (padj <", padj_cutoff, ", |log2FC| >=", log2FC, "):", length(unique(res_df_all_sig$gene)), "\n")

  if (run_both_tests) {
    lrt_sig <- nrow(res_df_all_sig[res_df_all_sig$DESeq_test == "LRT", ])
    wald_sig <- nrow(res_df_all_sig[res_df_all_sig$DESeq_test == "Wald", ])
    cat("  - LRT significant:", lrt_sig, "\n")
    cat("  - Wald significant:", wald_sig, "\n")
  }

  if (nrow(res_df_all_sig) > 0) {
    direction_counts <- table(res_df_all_sig$direction)
    cat("  - Upregulated:", ifelse("Up" %in% names(direction_counts), direction_counts["Up"], 0), "\n")
    cat("  - Downregulated:", ifelse("Down" %in% names(direction_counts), direction_counts["Down"], 0), "\n")
  }

  return(list(
    all_results = res_df_all,
    sig_results = res_df_all_sig,
    dds_wald = if (exists("dds_Wald")) dds_Wald else NULL,
    dds_lrt = if (run_both_tests && exists("dds_LRT")) dds_LRT else NULL,
    celltype = celltype,
    comparison = paste(diagnosis_group, "vs", reference_group),
    n_total = length(unique(res_df_all$gene)),
    n_significant = length(unique(res_df_all_sig$gene)),
    files = list(
      all_results = all_results_file,
      sig_results = sig_results_file
    )
  ))
}
