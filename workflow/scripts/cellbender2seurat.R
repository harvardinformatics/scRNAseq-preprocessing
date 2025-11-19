args <- commandArgs(trailingOnly = TRUE)
cellbender_h5 <- args[1]
output <- args[2]

library("Seurat")

h5_loaded <- Read10X_h5(cellbender_h5)
seurat <- CreateSeuratObject(counts = h5_loaded)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "(?i)^mt-")
seurat <- SCTransform(seurat, vars.to.regress = "percent.mt", verbose = FALSE)
seurat <- RunPCA(seurat, verbose = FALSE)
seurat <- RunUMAP(seurat, dims = 1:30)
seurat <- FindNeighbors(seurat, dims = 1:30)
seurat <- FindClusters(seurat)
saveRDS(seurat,file=output)

