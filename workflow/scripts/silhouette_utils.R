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

  message(sprintf(
    "Calculating approximate silhouette widths for %d cells with bluster::approxSilhouette",
    sum(valid_cells)
  ))

  sil <- bluster::approxSilhouette(
    emb[valid_cells, dims, drop = FALSE],
    clust[valid_cells]
  )
  sil_values[valid_cells] <- as.numeric(sil$width)
  seurat_obj@meta.data[[sil_col_name]] <- sil_values

  seurat_obj
}
