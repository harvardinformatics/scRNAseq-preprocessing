args <- commandArgs(trailingOnly = TRUE)
seurat <- args[1]
filtered_output <- args[2]
min_nfeature <- as.numeric(args[3])
min_ncount <- as.numeric(args[4])
max_mtdna <- as.numeric(args[5])
markers_output <- args[6]

library("Seurat")
library("glmGamPoi")
library("scater")
library("tidyverse")
library("tools")

options(future.globals.maxSize = 16 * 1024^3)

seurat <- readRDS(seurat)

seurat_filtered <- subset(seurat, subset = nFeature_RNA > min_nfeature & nCount_RNA > min_ncount & percent.mt < max_mtdna)

seurat_filtered <- SCTransform(seurat_filtered, vars.to.regress = "percent.mt", verbose = FALSE)
seurat_filtered <- RunPCA(seurat_filtered, verbose = FALSE)
seurat_filtered <- RunUMAP(seurat_filtered, dims = 1:30)
seurat_filtered <- FindNeighbors(seurat_filtered, dims = 1:30)
seurat_filtered <- FindClusters(seurat_filtered)
saveRDS(seurat_filtered,file=filtered_output)

markers <- FindAllMarkers(seurat_filtered)
markers$genesymbol <- row.names(markers)
workflow_name <- file_path_sans_ext(basename(markers_output))
workflow_name <- gsub("_markergenes", "", workflow_name)
sig_markers <- as_tibble(markers) %>% filter(p_val_adj<=0.05) %>% mutate(workflow=workflow_name)
write_csv(sig_markers,file=markers_output)
