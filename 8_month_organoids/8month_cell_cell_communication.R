#!/usr/bin/env Rscript

# Cell-cell communication using CellChat
# Open the conda environment with the "CellChat" package

library(Seurat)
library(CellChat)
library(dplyr)

astro_mo8 <- readRDS("/raidixshare_logg01/jchen/project/01_sarah_proj/03_seurat_files/seurat_mo8/seurat_8mo_annots_0411.rds")
DimPlot(astro_mo8, reduction = "umap", label = T, group.by = "cell_type") #just check
# Prepare subset - only focus on the cell types that we are interested in
unique(astro_mo8@meta.data$cell_type)
ccc_mo8 <- subset(astro_mo8, subset = cell_type %in% c("ASC", "OPC_ODC", "GABA_INH", "NKCC1_pos_cells", "KCC2_EX"))
# check distribution of different cell types
table(ccc_mo8@meta.data$patient_id, ccc_mo8@meta.data$cell_type)
table(ccc_mo8@meta.data$diagnosis, ccc_mo8@meta.data$cell_type)
# ccc
Obj1 = subset(ccc_mo8, subset = diagnosis == 'AD')
Obj2 = subset(ccc_mo8, subset = diagnosis == 'NC')
Matrix1 = Obj1@assays$RNA$data %>% na.omit 
Matrix2 = Obj2@assays$RNA$data %>% na.omit
cellchat_1 <- createCellChat(object = Matrix1, meta = Obj1@meta.data, group.by = "cell_type")
cellchat_2 <- createCellChat(object = Matrix2, meta = Obj2@meta.data, group.by = "cell_type")
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

### 1 - cellchat_1 - AD
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
# Other options:
# Only uses the Secreted Signaling from CellChatDB v1
# CellChatDB.use <- subsetDB(CellChatDB, search = list(c("Secreted Signaling"), c("CellChatDB v1")), key = c("annotation", "version"))

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB)

# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 

# set the used database in the object
cellchat_1@DB <- CellChatDB.use

ptm = Sys.time()
# subset the expression data of signaling genes for saving computation cost
cellchat_1 <- subsetData(cellchat_1) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
cellchat_1 <- identifyOverExpressedGenes(cellchat_1)
cellchat_1 <- identifyOverExpressedInteractions(cellchat_1)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 692

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
#> [1] 210.2585
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

ptm = Sys.time()
cellchat_1 <- computeCommunProb(cellchat_1, type = "triMean")
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2025-04-30 23:09:31.186255]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2025-04-30 23:16:10.464598]"
cellchat_1 <- filterCommunication(cellchat_1, min.cells = 10)
cellchat_1 <- computeCommunProbPathway(cellchat_1)
cellchat_1 <- aggregateNet(cellchat_1)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
#> [1] 466.5054

groupSize_1 <- as.numeric(table(cellchat_1@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_1@net$count, vertex.weight = groupSize_1, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_1@net$weight, vertex.weight = groupSize_1, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat_1 <- cellchat_1@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat_1)) {
  mat2 <- matrix(0, nrow = nrow(mat_1), ncol = ncol(mat_1), dimnames = dimnames(mat_1))
  mat2[i, ] <- mat_1[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize_1, weight.scale = T, edge.weight.max = max(mat_1), title.name = rownames(mat_1)[i])
}

cellchat_1@netP$pathways
#> 'MK''PTN''BMP''KIT''EGF'

# 'MK' pathway
pathways.show <- c("MK") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat_1, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_1, signaling = pathways.show, layout = "circle")

# 'BMP' pathway
pathways.show <- c("BMP")
netVisual_aggregate(cellchat_1, signaling = pathways.show,  vertex.receiver = vertex.receiver)
par(mfrow=c(1,1)) # Circle plot
netVisual_aggregate(cellchat_1, signaling = pathways.show, layout = "circle")

# 'PTN' pathway
pathways.show <- c("PTN") 
netVisual_aggregate(cellchat_1, signaling = pathways.show,  vertex.receiver = vertex.receiver)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_1, signaling = pathways.show, layout = "circle")

### 2 - cellchat_2 - NC
cellchat_2@DB <- CellChatDB.use
cellchat_2 <- subsetData(cellchat_2) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat_2 <- identifyOverExpressedGenes(cellchat_2)
cellchat_2 <- identifyOverExpressedInteractions(cellchat_2)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
ptm = Sys.time()
cellchat_2 <- computeCommunProb(cellchat_2, type = "triMean")
cellchat_2 <- filterCommunication(cellchat_2, min.cells = 10)
cellchat_2 <- computeCommunProbPathway(cellchat_2)
cellchat_2 <- aggregateNet(cellchat_2)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
cellchat_2 <- aggregateNet(cellchat_2)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
#> [1] 1023.098
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2025-04-30 23:35:01.361573]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2025-04-30 23:41:04.804629]"
#> [1] 373.449
#> [1] 373.4752

groupSize_2 <- as.numeric(table(cellchat_2@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_2@net$count, vertex.weight = groupSize_2, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_2@net$weight, vertex.weight = groupSize_2, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

cellchat_2@netP$pathways
#> 'PTN''MK''NRG''BMP''GRN'

# 'MK' pathway
pathways.show <- c("MK") 
netVisual_aggregate(cellchat_2, signaling = pathways.show,  vertex.receiver = vertex.receiver)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_2, signaling = pathways.show, layout = "circle")

# 'BMP' pathway
pathways.show <- c("BMP") 
netVisual_aggregate(cellchat_2, signaling = pathways.show,  vertex.receiver = vertex.receiver)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_2, signaling = pathways.show, layout = "circle")

# 'PTN' pathway
pathways.show <- c("PTN") 
netVisual_aggregate(cellchat_2, signaling = pathways.show,  vertex.receiver = vertex.receiver)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_2, signaling = pathways.show, layout = "circle")
