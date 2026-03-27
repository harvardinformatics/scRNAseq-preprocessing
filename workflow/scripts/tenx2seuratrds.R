args <- commandArgs(trailingOnly = TRUE)
filtered <- args[1]
output <- args[2]
markers_output<- args[3]

library("Seurat")
library("tidyverse")
library("glmGamPoi")
library("tools")

options(future.globals.maxSize = 16 * 1024^3)

filtered_loaded <- Seurat::Read10X(filtered)
seurat_obj <- CreateSeuratObject(counts = filtered_loaded)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "(?i)^mt-")
seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj)
saveRDS(seurat_obj,file=output)

markers <- FindAllMarkers(seurat_obj)
markers$genesymbol <- row.names(markers)
workflow_name <- file_path_sans_ext(basename(markers_output))
workflow_name <- gsub("_markergenes", "", workflow_name)
sig_markers <- as_tibble(markers) %>% filter(p_val_adj<=0.05) %>% mutate(workflow=workflow_name)
write_csv(sig_markers,file=markers_output)
