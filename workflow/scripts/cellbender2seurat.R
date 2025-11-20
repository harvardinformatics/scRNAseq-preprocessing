args <- commandArgs(trailingOnly = TRUE)
cellbender_h5 <- args[1]
output <- args[2]

library("Seurat")
library("scCustomize")
mat <- Read_CellBender_h5_Mat(cellbender_h5)
seurat <- CreateSeuratObject(mat)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "(?i)^mt-")
seurat <- SCTransform(seurat, vars.to.regress = "percent.mt", verbose = FALSE)
seurat <- RunPCA(seurat, verbose = FALSE)
seurat <- RunUMAP(seurat, dims = 1:30)
seurat <- FindNeighbors(seurat, dims = 1:30)
seurat <- FindClusters(seurat)
saveRDS(seurat,file=output)

