args <- commandArgs(trailingOnly = TRUE)
filtered <- args[1]
output <- args[2]
library("Seurat")
library("glmGamPoi")
filtered_loaded <- Seurat::Read10X(filtered)
seurat_obj <- CreateSeuratObject(counts = filtered_loaded)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "(?i)^mt-")
seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj)
saveRDS(seurat_obj,file=output)
