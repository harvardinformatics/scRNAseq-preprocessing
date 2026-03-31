args <- commandArgs(trailingOnly = TRUE)
cellbender_h5 <- args[1]
output <- args[2]
markers_output <- args[3]

library("Seurat")
library("tidyverse")
library("scCustomize")
library(tools)
options(future.globals.maxSize = 16 * 1024^3)

mat <- Read_CellBender_h5_Mat(cellbender_h5)
seurat <- CreateSeuratObject(mat)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "(?i)^mt-")
seurat <- SCTransform(seurat, vars.to.regress = "percent.mt", verbose = FALSE)
seurat <- RunPCA(seurat, verbose = FALSE)
seurat <- RunUMAP(seurat, dims = 1:30)
seurat <- FindNeighbors(seurat, dims = 1:30)
seurat <- FindClusters(seurat)
saveRDS(seurat,file=output)

markers <- FindAllMarkers(seurat)
markers$genesymbol <- row.names(markers)
workflow_name <- file_path_sans_ext(basename(markers_output))
workflow_name <- gsub("_markergenes", "", workflow_name)
sig_markers <- as_tibble(markers) %>% filter(p_val_adj<=0.05) %>% mutate(workflow=workflow_name)
write_csv(sig_markers,file=markers_output)
