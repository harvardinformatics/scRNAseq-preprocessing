args <- commandArgs(trailingOnly = TRUE)
seurat <- args[1]
filtered_output <- args[2]
markers_output <- args[3]

library("Seurat")
library("glmGamPoi")
library("scater")
library("tidyverse")
library("tools")

options(future.globals.maxSize = 16 * 1024^3)

seurat <- readRDS(seurat)
sce <-as.SingleCellExperiment(seurat)
mito_genes <- grep("(?i)^mt-", rownames(sce), value = TRUE)
qc <- perCellQCMetrics(sce, subsets = list(Mito = mito_genes))
high_mito <- isOutlier(qc$subsets_Mito_percent, nmads=3, type="higher")
low_umi     <- isOutlier(qc$sum, nmads = 3, type = "lower")
low_feature <- isOutlier(qc$detected, nmads = 3, type = "lower")
discard <- high_mito | low_umi | low_feature
cells_to_keep <- colnames(seurat)[!discard]
seurat_filtered<- subset(seurat, cells = cells_to_keep)

seurat_filtered[["percent.mt"]] <- PercentageFeatureSet(seurat_filtered, pattern = "(?i)^mt-")
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
