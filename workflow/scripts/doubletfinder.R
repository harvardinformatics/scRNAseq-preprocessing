args <- commandArgs(trailingOnly = TRUE)
seurat <- args[1]
output <- args[2]
nclusters_output <- args[3]
cluster_ids_output <- args[4]

library("Seurat")
library("DoubletFinder")
library("igraph")
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

sweep.res.list <- paramSweep(seurat, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
optimal_pk <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric),]$pK))

# set value of nExp
nExp_poi <- round(0.15*nrow(seurat@meta.data))
homotypic.prop <- modelHomotypic(seurat$seurat_clusters)
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

seurat <- doubletFinder(seurat, PCs = 1:10, pN = 0.25,
                                 pK = optimal_pk, nExp = nExp_poi.adj,
                                 reuse.pANN = NULL,
                                 sct=TRUE)

classifier_colname <- paste("DF.classifications_0.25",optimal_pk,nExp_poi.adj,sep="_")
seurat_singlets <- subset(seurat, cells = rownames(seurat@meta.data)[seurat@meta.data[[classifier_colname]] == "Singlet"]) 

seurat_singlets[["percent.mt"]] <- PercentageFeatureSet(seurat_singlets, pattern = "(?i)^mt-")
seurat_singlets <- SCTransform(seurat_singlets, vars.to.regress = "percent.mt", verbose = FALSE)
seurat_singlets <- RunPCA(seurat_singlets, verbose = FALSE)
seurat_singlets <- RunUMAP(seurat_singlets, dims = 1:30)
seurat_singlets <- FindNeighbors(seurat_singlets, dims = 1:30)
seurat_singlets <- FindClusters(seurat_singlets)
seurat_singlets <- add_silhouette_to_metadata(seurat_singlets)

saveRDS(seurat_singlets,file=output)
write_cluster_metadata(seurat_singlets, nclusters_output, cluster_ids_output)
