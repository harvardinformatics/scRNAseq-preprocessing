source("workflow/scripts/regression_utils.R")
opts <- parse_script_args(commandArgs(trailingOnly = TRUE))
args <- opts$positional
if (length(args) < 3) {
  stop("Usage: downsample_clusters.R <seurat_rds> <output_tsv> <n_replicates> [downsample_rate] ",
       "[--vars-to-regress=A,B] [--sample=NAME]", call. = FALSE)
}

seurat_rds <- args[1]
output <- args[2]
n_replicates <- suppressWarnings(as.integer(args[3]))
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
options(future.globals.maxSize = 16 * 1024^3)

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

SubSampleReSCTSeuratObject <- function(seurat_obj, subrate, replicate_seed, vars_to_regress) {
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

  subsampled_seurat <- SCTransform(subsampled_seurat, vars.to.regress = vars_to_regress, verbose = FALSE)
  subsampled_seurat <- RunPCA(subsampled_seurat, verbose = FALSE)
  pca_dims <- seq_len(min(30, ncol(Embeddings(subsampled_seurat, "pca"))))
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

  dat <- tibble::tibble(
    cell_id = names(seurat_obj@active.ident),
    cluster = seurat_obj$seurat_clusters
  ) %>%
    tidyr::nest(data = -cluster) %>%
    dplyr::arrange(cluster)

  for (original_cluster in unique(seurat_obj$presub_clusters)) {
    barcodes <- rownames(
      subset(seurat_obj@meta.data, presub_clusters == original_cluster)
    )

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

# Each replicate re-runs SCTransform/PCA/clustering from raw counts and only needs the
# original cluster labels for the Jaccard comparison. The input also carries a full SCT
# assay (a dense scale.data), PCA/UMAP embeddings and neighbor graphs from upstream
# clustering; those would be pinned for the whole loop and re-copied into every subset(),
# yet are never used here. Strip to a minimal counts-only object (measured ~36 -> ~30 GB
# peak per replicate on a 77k-cell dataset; the gap grows with cell count). Results are
# unchanged: SCTransform operates on the RNA counts, which are preserved exactly.
# Each replicate's SCTransform regresses out --vars-to-regress (none by default; see regression_utils.R).
# nCount_RNA and nFeature_RNA are recomputed below. percent.mt is the input's, which no step here uses unless it
# is regressed out; then it is recomputed on these RNA counts, as add_percent_mt() computes it in every step,
# from the mitochondrial features the input recorded from the sample's mt_genes list or those matching
# "(?i)^mt-" (a cell's value does not change when cells are subsampled).
mt_features <- stored_mt_features(seurat_obj)
keep_meta <- intersect(c("seurat_clusters", "percent.mt"), colnames(seurat_obj@meta.data))
seurat_obj <- CreateSeuratObject(
  counts = GetAssayData(seurat_obj, assay = "RNA", layer = "counts"),
  meta.data = seurat_obj@meta.data[, keep_meta, drop = FALSE]
)
gc(verbose = FALSE)
if ("percent.mt" %in% opts$vars_to_regress) {
  seurat_obj <- add_percent_mt(seurat_obj, mt_features, opts, context = "downsample_clusters")
}
check_regression_vars(seurat_obj, opts, context = "downsample_clusters")

replicate_results <- vector("list", n_replicates)
for (replicate in seq_len(n_replicates)) {
  t0 <- Sys.time()
  replicate_seed <- workflow_seed + replicate
  subsampled_obj <- SubSampleReSCTSeuratObject(seurat_obj, downsample_rate, replicate_seed, opts$vars_to_regress)
  replicate_results[[replicate]] <- GetJaccardMaxByCluster(subsampled_obj, replicate)
  rm(subsampled_obj)
  gc(verbose = FALSE)
  message(sprintf(
    "[%s] replicate %d/%d done: %.1f sec",
    format(Sys.time(), "%H:%M:%S"), replicate, n_replicates,
    as.numeric(Sys.time() - t0, units = "secs")
  ))
}

jaccard_max_stats <- dplyr::bind_rows(replicate_results)

dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)
write_tsv(jaccard_max_stats, output)
