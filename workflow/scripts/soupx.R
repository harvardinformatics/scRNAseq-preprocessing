args <- commandArgs(trailingOnly = TRUE)
filtered <- args[1]
raw <- args[2]
seurat_base <- args[3]
output <- args[4]
markers_output <- args[5]

library("Seurat")
library("glmGamPoi")
library("SoupX")
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

filtered_matrix <- Seurat::Read10X(filtered)
raw_matrix <- Seurat::Read10X(raw)


seurat_base <- readRDS(seurat_base)
soup_channel <- SoupX::SoupChannel(tod = raw_matrix,toc=filtered_matrix,
                is10X = TRUE)
soup_channel$tod <- raw
soup_channel <- SoupX::setClusters(soup_channel, 
                                   clusters = as.factor(Idents(seurat_base)))
soup_channel <- setDR(soup_channel, 
                DR=Seurat::Embeddings(seurat_base, "umap"))
soup_channel <- autoEstCont(soup_channel)
corrected_counts <- adjustCounts(soup_channel,roundToInt=TRUE)
seurat_soupx <- CreateSeuratObject(counts = corrected_counts)
seurat_soupx[["percent.mt"]] <- PercentageFeatureSet(seurat_soupx, pattern = "(?i)^mt-")
seurat_soupx <- SCTransform(seurat_soupx, vars.to.regress = "percent.mt", verbose = FALSE)
seurat_soupx <- RunPCA(seurat_soupx, verbose = FALSE)
seurat_soupx <- RunUMAP(seurat_soupx, dims = 1:30)
seurat_soupx <- FindNeighbors(seurat_soupx, dims = 1:30)
seurat_soupx <- FindClusters(seurat_soupx)
seurat_soupx <- add_silhouette_to_metadata(seurat_soupx)
saveRDS(seurat_soupx,file=output)

markers <- FindAllMarkers(seurat_soupx)
markers$genesymbol <- row.names(markers)
workflow_name <- file_path_sans_ext(basename(markers_output))
workflow_name <- gsub("_markergenes", "", workflow_name)
sig_markers <- as_tibble(markers) %>% filter(p_val_adj<=0.05) %>% mutate(workflow=workflow_name)
write_csv(sig_markers,file=markers_output)
