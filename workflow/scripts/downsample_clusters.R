args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: downsample_clusters.R <seurat_rds> <output_tsv> [n_replicates] [downsample_rate]", call. = FALSE)
}

seurat_rds <- args[1]
output <- args[2]
n_replicates <- if (length(args) >= 3) as.integer(args[3]) else as.integer(Sys.getenv("SCRNASEQ_DOWNSAMPLE_REPLICATES", "100"))
downsample_rate <- if (length(args) >= 4) as.numeric(args[4]) else as.numeric(Sys.getenv("SCRNASEQ_DOWNSAMPLE_RATE", "0.8"))
workflow_seed <- suppressWarnings(as.integer(Sys.getenv("SCRNASEQ_DOWNSAMPLE_SEED", "12345")))

if (length(n_replicates) != 1 || is.na(n_replicates) || n_replicates < 1) {
  stop("n_replicates must be a positive integer", call. = FALSE)
}
if (length(downsample_rate) != 1 || is.na(downsample_rate) || downsample_rate <= 0 || downsample_rate > 1) {
  stop("downsample_rate must be > 0 and <= 1", call. = FALSE)
}
if (length(workflow_seed) != 1 || is.na(workflow_seed)) {
  stop("SCRNASEQ_DOWNSAMPLE_SEED must be an integer", call. = FALSE)
}
set.seed(workflow_seed)

suppressPackageStartupMessages({
  library("tidyverse")
  library("Seurat")
  library("glmGamPoi")
})
options(future.globals.maxSize = 2 * 1024^3)

JaccardSimilarity <- function(set1, set2) {
  intersect_length <- length(intersect(set1, set2))
  union_length <- length(set1) + length(set2) - intersect_length
  intersect_length / union_length
}

RandomSubsetData <- function(object, rate, random.subset.seed = NULL, ...) {
  ncells <- nrow(object@meta.data)
  ncells.subsample <- round(ncells * rate)

  set.seed(random.subset.seed)

  selected.cells <- sample(colnames(object), ncells.subsample)
  object <- subset(object, cells = selected.cells, ...)
  return(object)
}

SubSampleReSCTSeuratObject <- function(seurat_obj, subrate, replicate_seed) {
  set.seed(replicate_seed)
  subsampled_seurat <- RandomSubsetData(
    seurat_obj,
    rate = subrate,
    random.subset.seed = replicate_seed
  )
  subsampled_seurat$presub_clusters <- seurat_obj@meta.data[
    colnames(subsampled_seurat),
    "seurat_clusters",
    drop = TRUE
  ]

  sct_args <- list(object = subsampled_seurat, verbose = FALSE)
  if ("percent.mt" %in% colnames(subsampled_seurat@meta.data)) {
    sct_args$vars.to.regress <- "percent.mt"
  }
  subsampled_seurat <- do.call(SCTransform, sct_args)
  subsampled_seurat <- RunPCA(subsampled_seurat, verbose = FALSE)
  pca_dims <- seq_len(min(30, ncol(Embeddings(subsampled_seurat, "pca"))))
  subsampled_seurat <- RunUMAP(
    subsampled_seurat,
    dims = pca_dims,
    seed.use = replicate_seed,
    verbose = FALSE
  )
  subsampled_seurat <- FindNeighbors(subsampled_seurat, dims = pca_dims, verbose = FALSE)
  subsampled_seurat <- FindClusters(
    subsampled_seurat,
    random.seed = replicate_seed,
    verbose = FALSE
  )
  return(subsampled_seurat)
}

GetJaccardMaxByCluster <- function(seurat_obj, bootstrap) {
  jaccard_max_stats <- tibble::tibble(
    clusterid = factor(),
    max_jaccard = numeric(),
    bootstrap_number = integer()
  )

  for (original_cluster in unique(seurat_obj$presub_clusters)) {
    barcodes <- rownames(
      subset(seurat_obj@meta.data, presub_clusters == original_cluster)
    )

    dat <- tibble::tibble(
      cell_id = names(seurat_obj@active.ident),
      cluster = seurat_obj$seurat_clusters
    ) %>%
      tidyr::nest(data = -cluster) %>%
      dplyr::arrange(cluster)

    maxstat <- dat %>%
      dplyr::mutate(
        jaccard = purrr::map(data, ~ JaccardSimilarity(barcodes, .x$cell_id))
      ) %>%
      dplyr::pull(jaccard) %>%
      unlist() %>%
      max()

    jaccard_max_stats <- jaccard_max_stats %>%
      tibble::add_row(
        clusterid = original_cluster,
        max_jaccard = maxstat,
        bootstrap_number = bootstrap
      )
  }

  return(jaccard_max_stats)
}

seurat_obj <- readRDS(seurat_rds)
if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
  stop("Input Seurat object is missing required metadata column: seurat_clusters", call. = FALSE)
}

jaccard_max_stats <- tibble::tibble(
  clusterid = factor(),
  max_jaccard = numeric(),
  bootstrap_number = integer()
)

for (replicate in seq_len(n_replicates)) {
  replicate_seed <- workflow_seed + replicate
  subsampled_obj <- SubSampleReSCTSeuratObject(seurat_obj, downsample_rate, replicate_seed)
  replicate_stats <- GetJaccardMaxByCluster(subsampled_obj, replicate)
  jaccard_max_stats <- bind_rows(jaccard_max_stats, replicate_stats)
}

dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)
write_tsv(jaccard_max_stats, output)
