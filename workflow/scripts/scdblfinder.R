args <- commandArgs(trailingOnly = TRUE)
seurat <- args[1]
output <- args[2]
markers_output <- args[3]

library("Seurat")
library("scDblFinder")
library("tidyverse")
library("tools")

options(future.globals.maxSize = 16 * 1024^3)

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

markers <- FindAllMarkers(seurat_singlets)
markers$genesymbol <- row.names(markers)
workflow_name <- file_path_sans_ext(basename(markers_output))
workflow_name <- gsub("_markergenes", "", workflow_name)
sig_markers <- as_tibble(markers) %>% filter(p_val_adj<=0.05) %>% mutate(workflow=workflow_name)
write_csv(sig_markers,file=markers_output)
