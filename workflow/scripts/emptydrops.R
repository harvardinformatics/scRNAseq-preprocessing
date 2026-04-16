args <- commandArgs(trailingOnly = TRUE)
raw <- args[1]
seurat_output <- args[2]
matrix_outdir <- args[3]
nclusters_output <- args[4]
cluster_ids_output <- args[5]

library("Seurat")
library("DropletUtils")
library("scater")
library("glmGamPoi")

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

# for back conversion to matrix
library("Matrix")
library("R.utils")

# run emptyDrops
droputil_rawdata <-read10xCounts(raw)
colnames(droputil_rawdata) <- droputil_rawdata$Barcode
rownames(droputil_rawdata) <- make.unique(rowData(droputil_rawdata)$Symbol)
dropletutils_out <-emptyDrops(counts(droputil_rawdata))
is.cell <-!is.na(dropletutils_out$FDR) & dropletutils_out$FDR<0.01
cell_barcodes <- colnames(droputil_rawdata)[which(is.cell)]
droputil_filtered <- droputil_rawdata[,is.cell]
droputil_filtered <- logNormCounts(droputil_filtered)
seurat_droputil_filtered <- as.Seurat(droputil_filtered, counts = "counts",data="logcounts")
seurat_droputil_filtered <- RenameAssays(seurat_droputil_filtered,originalexp="RNA")
seurat_droputil_filtered[["percent.mt"]] <- PercentageFeatureSet(seurat_droputil_filtered, pattern = "(?i)^mt-")
seurat_droputil_filtered <- SCTransform(seurat_droputil_filtered, vars.to.regress = "percent.mt", verbose = FALSE)
seurat_droputil_filtered <- RunPCA(seurat_droputil_filtered, verbose = FALSE)
seurat_droputil_filtered <- RunUMAP(seurat_droputil_filtered, dims = 1:30)
seurat_droputil_filtered <- FindNeighbors(seurat_droputil_filtered, dims = 1:30)
seurat_droputil_filtered <- FindClusters(seurat_droputil_filtered)
seurat_droputil_filtered <- add_silhouette_to_metadata(seurat_droputil_filtered)
saveRDS(seurat_droputil_filtered,file=seurat_output)
write_cluster_metadata(seurat_droputil_filtered, nclusters_output, cluster_ids_output)


# convert emptyDrops-filtered object to 10x matrix data structure
cm <- counts(droputil_filtered)
if (!inherits(cm, "dgCMatrix")) {
    cm <- as(cm, "dgCMatrix")
}


dir.create(matrix_outdir, showWarnings = FALSE)

# 1. Write uncompressed MatrixMarket file
writeMM(cm, file.path(matrix_outdir, "matrix.mtx"))

# 2. Write barcodes
writeLines(colnames(droputil_filtered), file.path(matrix_outdir, "barcodes.tsv"))

# 3. Write features
features <- rownames(droputil_filtered)
feature_df <- data.frame(
    id   = features,
    name = features,
    feature_type = "Gene Expression"
)
write.table(
    feature_df,
    file.path(matrix_outdir, "features.tsv"),
    sep = "\t",
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE
)

# 4. Gzip them using R's built-in compression
R.utils::gzip(file.path(matrix_outdir, "matrix.mtx"), overwrite = TRUE)
R.utils::gzip(file.path(matrix_outdir, "barcodes.tsv"), overwrite = TRUE)
R.utils::gzip(file.path(matrix_outdir, "features.tsv"), overwrite = TRUE)
