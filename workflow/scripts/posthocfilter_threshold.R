args <- commandArgs(trailingOnly = TRUE)
seurat <- args[1]
filtered_output <- args[2]
min_nfeature <- as.numeric(args[3])
min_ncount <- as.numeric(args[4])
max_mtdna <- as.numeric(args[5])
nclusters_output <- args[6]
cluster_ids_output <- args[7]

library("Seurat")
library("glmGamPoi")
library("scater")

options(future.globals.maxSize = 16 * 1024^3)

source("workflow/scripts/silhouette_utils.R")

write_cluster_metadata <- function(seurat_obj, nclusters_output, cluster_ids_output) {
  cluster_ids <- levels(Idents(seurat_obj))
  if (is.null(cluster_ids) || length(cluster_ids) == 0) {
    cluster_ids <- sort(unique(as.character(Idents(seurat_obj))))
  }
  cluster_ids <- cluster_ids[!is.na(cluster_ids) & nzchar(cluster_ids)]
  writeLines(cluster_ids, con = cluster_ids_output)
  writeLines(as.character(length(cluster_ids)), con = nclusters_output)
}

seurat <- readRDS(seurat)

seurat_filtered <- subset(seurat, subset = nFeature_RNA > min_nfeature & nCount_RNA > min_ncount & percent.mt < max_mtdna)

seurat_filtered <- SCTransform(seurat_filtered, vars.to.regress = "percent.mt", verbose = FALSE)
seurat_filtered <- RunPCA(seurat_filtered, verbose = FALSE)
seurat_filtered <- RunUMAP(seurat_filtered, dims = 1:30)
seurat_filtered <- FindNeighbors(seurat_filtered, dims = 1:30)
seurat_filtered <- FindClusters(seurat_filtered)
seurat_filtered <- add_silhouette_to_metadata(seurat_filtered)
saveRDS(seurat_filtered,file=filtered_output)
write_cluster_metadata(seurat_filtered, nclusters_output, cluster_ids_output)
