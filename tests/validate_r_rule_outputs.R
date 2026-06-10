args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("usage: validate_r_rule_outputs.R SEURAT_RDS MARKER_CSV", call. = FALSE)
}

rds_path <- args[[1]]
marker_path <- args[[2]]

suppressPackageStartupMessages({
    library(SeuratObject)
})

fail <- function(message) {
    stop(message, call. = FALSE)
}

if (!file.exists(rds_path) || file.info(rds_path)$size == 0) {
    fail(paste("missing or empty Seurat RDS:", rds_path))
}
if (!file.exists(marker_path) || file.info(marker_path)$size == 0) {
    fail(paste("missing or empty marker CSV:", marker_path))
}

obj <- readRDS(rds_path)
if (!inherits(obj, "Seurat")) {
    fail(paste("RDS is not a Seurat object; classes:", paste(class(obj), collapse = ",")))
}

n_features <- nrow(obj)
n_cells <- ncol(obj)
if (n_features < 100) {
    fail(paste("Seurat object has too few features:", n_features))
}
if (n_cells < 100) {
    fail(paste("Seurat object has too few cells:", n_cells))
}

metadata <- obj@meta.data
if (nrow(metadata) != n_cells) {
    fail("metadata row count does not match Seurat cell count")
}
if (anyDuplicated(rownames(metadata))) {
    fail("metadata contains duplicate cell barcodes")
}
if (any(is.na(rownames(metadata))) || any(rownames(metadata) == "")) {
    fail("metadata contains missing or empty cell barcodes")
}

required_metadata <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "seurat_clusters")
missing_metadata <- setdiff(required_metadata, colnames(metadata))
if (length(missing_metadata) > 0) {
    fail(paste("metadata missing columns:", paste(missing_metadata, collapse = ",")))
}

if (any(!is.finite(metadata$nCount_RNA)) || any(metadata$nCount_RNA <= 0)) {
    fail("metadata nCount_RNA contains non-positive or non-finite values")
}
if (any(!is.finite(metadata$nFeature_RNA)) || any(metadata$nFeature_RNA <= 0)) {
    fail("metadata nFeature_RNA contains non-positive or non-finite values")
}
if (any(!is.finite(metadata$percent.mt)) || any(metadata$percent.mt < 0 | metadata$percent.mt > 100)) {
    fail("metadata percent.mt contains values outside [0, 100]")
}

cluster_ids <- as.character(metadata$seurat_clusters)
if (length(unique(cluster_ids)) < 2) {
    fail("Seurat object has fewer than two clusters")
}

required_reductions <- c("pca", "umap")
missing_reductions <- setdiff(required_reductions, names(obj@reductions))
if (length(missing_reductions) > 0) {
    fail(paste("Seurat object missing reductions:", paste(missing_reductions, collapse = ",")))
}

markers <- read.csv(marker_path, stringsAsFactors = FALSE, check.names = FALSE)
expected_marker_columns <- c("genesymbol", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "workflow")
if (!identical(colnames(markers), expected_marker_columns)) {
    fail(paste(
        "marker CSV columns differ; observed:", paste(colnames(markers), collapse = ","),
        "expected:", paste(expected_marker_columns, collapse = ",")
    ))
}
if (nrow(markers) == 0) {
    fail("marker CSV has no rows")
}
if (any(is.na(markers$genesymbol)) || any(markers$genesymbol == "")) {
    fail("marker CSV contains missing or empty genesymbol values")
}

numeric_columns <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
for (column in numeric_columns) {
    values <- markers[[column]]
    if (any(!is.finite(values))) {
        fail(paste("marker CSV column has non-finite values:", column))
    }
}
for (column in c("p_val", "pct.1", "pct.2", "p_val_adj")) {
    values <- markers[[column]]
    if (any(values < 0 | values > 1)) {
        fail(paste("marker CSV column has values outside [0, 1]:", column))
    }
}

marker_clusters <- as.character(markers$cluster)
if (!all(marker_clusters %in% cluster_ids)) {
    fail("marker CSV contains clusters not present in Seurat metadata")
}
if (length(unique(marker_clusters)) < 2) {
    fail("marker CSV has markers for fewer than two clusters")
}
if (!all(markers$workflow == "filtered_seurat_tenx_test")) {
    fail("marker CSV workflow column does not match filtered_seurat_tenx_test")
}

cat(sprintf(
    "validated R rule outputs: %d features, %d cells, %d marker rows, %d clusters\n",
    n_features,
    n_cells,
    nrow(markers),
    length(unique(cluster_ids))
))
