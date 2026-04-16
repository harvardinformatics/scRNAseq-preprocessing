args <- commandArgs(trailingOnly = TRUE)
seurat <- args[1]
output <- args[2]
markers_output <- args[3]

library("Seurat")
library("scDblFinder")
library("tidyverse")
library("tools")

options(future.globals.maxSize = 16 * 1024^3)

add_silhouette_to_metadata <- function(
    seurat_obj,
    cluster_col = "seurat_clusters",
    reduction = "pca",
    dims = 1:30,
    sil_col_name = "silhouette_width"
) {
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", cluster_col, "not found in meta.data"))
  }
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(paste("Reduction", reduction, "not found in Seurat object"))
  }

  emb <- Embeddings(seurat_obj, reduction = reduction)
  dims <- dims[dims <= ncol(emb)]
  if (length(dims) == 0) {
    stop(paste("No", reduction, "dimensions available for silhouette calculation"))
  }

  clust <- seurat_obj@meta.data[[cluster_col]]
  valid_cells <- !is.na(clust)
  sil_values <- rep(NA_real_, length(clust))

  if (sum(valid_cells) < 2 || length(unique(clust[valid_cells])) < 2) {
    seurat_obj@meta.data[[sil_col_name]] <- sil_values
    return(seurat_obj)
  }

  clust_int <- as.integer(as.factor(clust[valid_cells]))
  sil <- cluster::silhouette(clust_int, stats::dist(emb[valid_cells, dims, drop = FALSE]))
  sil_values[valid_cells] <- sil[, "sil_width"]
  seurat_obj@meta.data[[sil_col_name]] <- sil_values

  seurat_obj
}

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
seurat_singlets <- add_silhouette_to_metadata(seurat_singlets)
saveRDS(seurat_singlets,file=output)

markers <- FindAllMarkers(seurat_singlets)
markers$genesymbol <- row.names(markers)
workflow_name <- file_path_sans_ext(basename(markers_output))
workflow_name <- gsub("_markergenes", "", workflow_name)
sig_markers <- as_tibble(markers) %>% filter(p_val_adj<=0.05) %>% mutate(workflow=workflow_name)
write_csv(sig_markers,file=markers_output)
