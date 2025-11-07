args <- commandArgs(trailingOnly = TRUE)
filtered <- args[1]

library("Seurat")

filtered_loaded <- Seurat::Read10X(filtered)
seurat <- CreateSeuratObject(counts = filtered_loaded)
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "(?i)^mt-")
seurat <- SCTransform(seurat, vars.to.regress = "percent.mt", verbose = FALSE)
seurat <- RunPCA(seurat, verbose = FALSE)
seurat <- RunUMAP(seurat, dims = 1:30)
seurat <- FindNeighbors(seurat, dims = 1:30)
seurat <- FindClusters(seurat)
saveRDS(seurat,file="
