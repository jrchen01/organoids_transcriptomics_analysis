# Unsupervised clustering
# Looping for different resolution

library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(ggrepel)

## step1 - Integration of the sequencing dataset across donor cell lines (patients)
astro.rna <- readRDS("/raidixshare_logg01/thais/SarahF/organoids/seurat_noMatrigel_0.05.rds")
raw.counts <- GetAssayData(astro.rna, layer = "counts")
raw.seurat <- CreateSeuratObject(counts = raw.counts)
metadata <- astro.rna@meta.data
raw.seurat <- AddMetaData(raw.seurat, metadata = metadata)
seurat.8mo <- subset(raw.seurat, subset = Time == "8 months") # extract 8 month data
seurat.8mo <- subset(seurat.8mo, subset = nFeature_RNA > 200)

seurat.8mo[["RNA"]] <- split(seurat.8mo[["RNA"]], f = seurat.8mo$patient_id)
seurat.8mo <- NormalizeData(seurat.8mo)
seurat.8mo <- FindVariableFeatures(seurat.8mo)
seurat.8mo <- ScaleData(seurat.8mo)
seurat.8mo <- RunPCA(seurat.8mo)

seurat.8mo <- FindNeighbors(seurat.8mo, dims = 1:30, reduction = "pca")
seurat.8mo <- FindClusters(seurat.8mo, resolution = 2, cluster.name = "unintegrated_clusters")

seurat.8mo <- RunUMAP(seurat.8mo, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(seurat.8mo, reduction = "umap.unintegrated", group.by = "patient_id") #unintegrated

seurat.8mo <- IntegrateLayers(object = seurat.8mo, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                              verbose = FALSE)
seurat.8mo[["RNA"]] <- JoinLayers(seurat.8mo[["RNA"]])
seurat.8mo #check, should be sth like "4 dimensional reductions calculated: pca, umap.unintegrated, integrated.cca, umap"

## step2 - loop for different resolution and different dims

seurat.8mo <- RunUMAP(seurat.8mo, dims = 1:20, reduction = "integrated.cca")
DimPlot(seurat.8mo, reduction = "umap", group.by = "patient_id") #initial check
ElbowPlot(seurat.8mo, reduction = "integrated.cca", ndims = 50) ##### needs to be confirmed

## step 2.1 - dims = 1:20
output_dir <- "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/plots/month8_seurat/diff_resol/umap/"
resolutions <- seq(0.05, 0.3, by = 0.05)
seurat.8mo <- FindNeighbors(seurat.8mo, reduction = "integrated.cca", dims = 1:20)

for (res in resolutions) {
  seurat.8mo <- FindClusters(seurat.8mo, resolution = res)
  seurat.8mo <- RunUMAP(seurat.8mo, dims = 1:20, reduction = "integrated.cca")
  p <- DimPlot(seurat.8mo, reduction = "umap") + ggtitle(paste("Resolution:", res, "| dim=20"))
  ggsave(filename = paste0(output_dir, "umap_dim20_resol_", res, ".png"), plot = p, width = 8, height = 6)
}

## step 2.2 - dims = 1:30
output_dir <- "/raidixshare_logg01/jchen/project/01_sarah_proj/01_jupyternotebook_files/plots/month8_seurat/diff_resol/umap/"
resolutions <- seq(0.05, 0.3, by = 0.05)
seurat.8mo <- FindNeighbors(seurat.8mo, reduction = "integrated.cca", dims = 1:30)

for (res in resolutions) {
  seurat.8mo <- FindClusters(seurat.8mo, resolution = res)
  seurat.8mo <- RunUMAP(seurat.8mo, dims = 1:30, reduction = "integrated.cca")
  p <- DimPlot(seurat.8mo, reduction = "umap") + ggtitle(paste("Resolution:", res, "| dim=30"))
  ggsave(filename = paste0(output_dir, "umap_dim30_resol_", res, ".png"), plot = p, width = 8, height = 6)
}

# choose dim = 30, reso = 0.15
astro.rna <- readRDS("/raidixshare_logg01/thais/SarahF/organoids/seurat_noMatrigel_0.05.rds")
raw.counts <- GetAssayData(astro.rna, layer = "counts")
raw.seurat <- CreateSeuratObject(counts = raw.counts)
metadata <- astro.rna@meta.data
raw.seurat <- AddMetaData(raw.seurat, metadata = metadata)
seurat.8mo <- subset(raw.seurat, subset = Time == "8 months")
seurat.8mo <- subset(seurat.8mo, subset = nFeature_RNA > 200)
seurat.8mo[["RNA"]] <- split(seurat.8mo[["RNA"]], f = seurat.8mo$patient_id)
seurat.8mo <- NormalizeData(seurat.8mo)
seurat.8mo <- FindVariableFeatures(seurat.8mo)
seurat.8mo <- ScaleData(seurat.8mo)
seurat.8mo <- RunPCA(seurat.8mo)
seurat.8mo <- FindNeighbors(seurat.8mo, dims = 1:30, reduction = "pca")
seurat.8mo <- FindClusters(seurat.8mo, resolution = 2, cluster.name = "unintegrated_clusters")
seurat.8mo <- RunUMAP(seurat.8mo, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
seurat.8mo <- IntegrateLayers(object = seurat.8mo, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                              verbose = FALSE)
seurat.8mo[["RNA"]] <- JoinLayers(seurat.8mo[["RNA"]])
seurat.8mo <- FindNeighbors(seurat.8mo, reduction = "integrated.cca", dims = 1:30)
seurat.8mo <- FindClusters(seurat.8mo, resolution = 0.15)
seurat.8mo <- RunUMAP(seurat.8mo, dims = 1:30, reduction = "integrated.cca")
DimPlot(seurat.8mo, reduction = "umap") # just check
output_dir <- "/raidixshare_logg01/jchen/project/01_sarah_proj/03_seurat_files/seurat_mo8/"
saveRDS(seurat.8mo, file = paste0(output_dir, "seurat_8mo_dim30_res0.15.rds"))

# check
astro_mo8 <- readRDS("/raidixshare_logg01/jchen/project/01_sarah_proj/03_seurat_files/seurat_mo8/seurat_8mo_dim30_res0.15.rds")
DimPlot(astro_mo8, reduction = "umap")
DimPlot(astro_mo8, reduction = "umap", group.by = "patient_id")
DimPlot(astro_mo8, reduction = "umap", group.by = "Phase")