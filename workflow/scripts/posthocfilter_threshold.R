args <- commandArgs(trailingOnly = TRUE)
seurat <- args[1]
filtered_output <- args[2]
min_nfeature <- args[3]
min_ncount <- args[4]
max_mtdna <- args[5]
library("Seurat")
library("glmGamPoi")
library("scater")

seurat <- readRDS(seurat)

seurat_filtered <- subset(seurat, subset = nFeature_RNA > min_nfeature & nCount_RNA > min_ncount & percent.mt < max_mtdna)

seurat_filtered <- SCTransform(seurat_filtered, vars.to.regress = "percent.mt", verbose = FALSE)
seurat_filtered <- RunPCA(seurat_filtered, verbose = FALSE)
seurat_filtered <- RunUMAP(seurat_filtered, dims = 1:30)
seurat_filtered <- FindNeighbors(seurat_filtered, dims = 1:30)
seurat_filtered <- FindClusters(seurat_filtered)
saveRDS(seurat_filtered,file=filtered_output)
