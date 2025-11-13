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
saveRDS(seurat_singlets,file=output)
