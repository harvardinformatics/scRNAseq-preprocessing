args <- commandArgs(trailingOnly = TRUE)
seurat <- args[1]
output <- args[2]

library("Seurat")
library("scDblFinder")

seurat <- readRDS(seurat)
sce <- as.SingleCellExperiment(seurat)
sce <- scDblFinder(sce)
seurat$scDblFinder.class <- colData(sce)$scDblFinder.class
seurat_singlets <- subset(seurat, subset = scDblFinder.class == "singlet")
seurat_singlets[["percent.mt"]] <- PercentageFeatureSet(seurat_singlets, pattern = "(?i)^mt-")
seurat_singlets <- SCTransform(seurat_singlets, vars.to.regress = "percent.mt", verbose = FALSE)
seurat_singlets <- RunPCA(seurat_singlets, verbose = FALSE)
seurat_singlets <- RunUMAP(seurat_singlets, dims = 1:30)
seurat_singlets <- FindNeighbors(seurat_singlets, dims = 1:30)
seurat_singlets <- FindClusters(seurat_singlets)
saveRDS(seurat_singlets,file=output)
