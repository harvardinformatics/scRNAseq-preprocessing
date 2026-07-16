# Generates a mid-scale downsampling test fixture by stratified resampling
# (with light Poisson jitter) of the existing tiny fixture, preserving
# per-cluster proportions. This exists purely to give the scaling test two
# data points of meaningfully different size, so quadratic-or-worse
# performance regressions (e.g. the do.call()/SCTransform hang) are
# detectable in CI without needing production-scale data.
#
# Usage: Rscript make_medium_fixture.R <input_rds> <output_rds> <target_ncells> [seed]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: make_medium_fixture.R <input_rds> <output_rds> <target_ncells> [seed]", call. = FALSE)
}

input_rds <- args[1]
output_rds <- args[2]
target_ncells <- as.integer(args[3])
seed <- if (length(args) >= 4) as.integer(args[4]) else 12345L

suppressPackageStartupMessages(library(Seurat))
set.seed(seed)

obj <- readRDS(input_rds)
counts <- obj[["RNA"]]$counts
clusters <- obj$seurat_clusters

cluster_props <- table(clusters) / length(clusters)
per_cluster_target <- round(cluster_props * target_ncells)

sampled_source_idx <- unlist(lapply(names(per_cluster_target), function(cl) {
  pool <- which(clusters == cl)
  sample(pool, size = per_cluster_target[[cl]], replace = TRUE)
}))

sampled_counts <- counts[, sampled_source_idx, drop = FALSE]
sampled_clusters <- unname(clusters[sampled_source_idx])

sampled_counts@x <- sampled_counts@x + rpois(length(sampled_counts@x), lambda = 0.02)

# Sparse jitter alone doesn't reliably de-duplicate resampled cells: a cell with few
# nonzero genes has a non-negligible chance every jitter draw is 0 (P(X=0) ~= 0.98
# per entry with lambda=0.02), leaving repeat draws of the same source cell exactly
# identical. Guarantee distinctness directly for any source index drawn more than
# once by bumping one entry per repeat draw.
repeat_draw <- which(duplicated(sampled_source_idx))
for (k in seq_along(repeat_draw)) {
  col <- repeat_draw[k]
  row <- ((k - 1) %% nrow(sampled_counts)) + 1
  sampled_counts[row, col] <- sampled_counts[row, col] + 1
}

colnames(sampled_counts) <- paste0("cell_", seq_len(ncol(sampled_counts)))

new_obj <- CreateSeuratObject(counts = sampled_counts)
pct_mt <- PercentageFeatureSet(new_obj, pattern = "^MT-")
new_obj$percent.mt <- if (is.data.frame(pct_mt)) pct_mt[[1]] else unname(pct_mt)
new_obj$seurat_clusters <- factor(sampled_clusters, levels = levels(clusters))
Idents(new_obj) <- new_obj$seurat_clusters

dir.create(dirname(output_rds), showWarnings = FALSE, recursive = TRUE)
saveRDS(new_obj, output_rds)

message(sprintf(
  "Wrote %s: %d cells x %d genes (%d clusters)",
  output_rds, ncol(new_obj), nrow(new_obj), length(levels(new_obj$seurat_clusters))
))
