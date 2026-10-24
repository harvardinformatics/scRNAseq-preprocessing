args <- commandArgs(trailingOnly = TRUE)
filtered <- args[1]
raw <- args[2]              # path to the raw 10x matrix dir (the background / empty-droplet source)
seurat_base <- args[3]
output <- args[4]
nclusters_output <- args[5]
cluster_ids_output <- args[6]

library("Seurat")
library("glmGamPoi")
library("decontX")

options(future.globals.maxSize = 16 * 1024^3)

# Cap the decontX EM's RcppParallel threads to the job's allocation. decontX exposes no
# thread argument, and RcppParallel otherwise defaults to every core visible on the node
# (not the cpus SLURM granted this task), which would oversubscribe a shared node. Harmless
# if the compiled core does not actually use RcppParallel. See DECONTX_THREADS in decontx.smk.
try(
  RcppParallel::setThreadOptions(numThreads = as.integer(Sys.getenv("DECONTX_THREADS", "1"))),
  silent = TRUE
)

source("workflow/scripts/silhouette_utils.R")
WORKFLOW_SEED <- set_workflow_seed()

write_cluster_metadata <- function(seurat_obj, nclusters_output, cluster_ids_output) {
  cluster_ids <- levels(Idents(seurat_obj))
  if (is.null(cluster_ids) || length(cluster_ids) == 0) {
    cluster_ids <- sort(unique(as.character(Idents(seurat_obj))))
  }
  cluster_ids <- cluster_ids[!is.na(cluster_ids) & nzchar(cluster_ids)]
  writeLines(cluster_ids, con = cluster_ids_output)
  writeLines(as.character(length(cluster_ids)), con = nclusters_output)
}

filtered_matrix <- Seurat::Read10X(filtered)
seurat_base <- readRDS(seurat_base)

# Reuse the base Seurat object's clusters as decontX's cell populations (z), exactly as the
# soupx arm reuses them for SoupX. decontX estimates each cell's contamination from the other
# populations' expression profiles, so a shared clustering keeps the arms comparable and avoids
# decontX's slower internal dbscan step. Align the clusters to the filtered matrix columns.
base_clusters <- as.factor(Idents(seurat_base))
common <- intersect(colnames(filtered_matrix), names(base_clusters))
if (length(common) == 0) {
  stop("decontx: no barcodes shared between the filtered matrix and the base Seurat object")
}
decon_counts <- filtered_matrix[, common, drop = FALSE]
z <- droplevels(base_clusters[common])
rm(seurat_base)
gc()

# decontX uses the raw matrix as the background (the SoupX-equivalent mode): the ambient profile
# is fixed to the empty-droplet gene totals. decontX does not align genes between x and the
# background, so the raw matrix is subset and reordered to the filtered matrix's genes; decontX
# itself drops any background barcode that also appears among the called cells, so the background
# resolves to the empty droplets.
raw_matrix <- Seurat::Read10X(raw)
background <- raw_matrix[rownames(decon_counts), , drop = FALSE]

# decontX requires >= 2 cell populations when z is supplied; fall back to its internal
# clustering for the rare single-cluster base object rather than erroring out.
z_use <- if (nlevels(z) >= 2) as.integer(z) else NULL

set.seed(WORKFLOW_SEED)
decontx_result <- decontX(
  x = decon_counts,
  z = z_use,
  background = background,
  seed = WORKFLOW_SEED,
  verbose = FALSE
)

# decontX returns fractional counts; round to integers for the downstream integer-count tools
# (SCTransform, the doublet callers) and drop the zeros that rounding introduces.
corrected_counts <- Matrix::drop0(round(decontx_result$decontXcounts))
# decontX + integer rounding can strip a low-count, high-contamination cell down to zero total
# counts; an all-zero column gives log_umi = log(0) = -Inf and aborts SCTransform. Drop any such
# cells (they carry no signal after decontamination) before building the Seurat object.
nonempty_cells <- Matrix::colSums(corrected_counts) > 0
if (!all(nonempty_cells)) {
  message(sprintf("decontx: dropping %d cell(s) reduced to zero counts after decontamination",
                  sum(!nonempty_cells)))
  corrected_counts <- corrected_counts[, nonempty_cells, drop = FALSE]
}
rm(filtered_matrix, decon_counts)
rm(raw_matrix, background)
rm(decontx_result)
gc()

seurat_decontx <- CreateSeuratObject(counts = corrected_counts)
rm(corrected_counts)
gc()
seurat_decontx[["percent.mt"]] <- PercentageFeatureSet(seurat_decontx, pattern = "(?i)^mt-")
seurat_decontx <- SCTransform(seurat_decontx, vars.to.regress = "percent.mt", seed.use = WORKFLOW_SEED, verbose = FALSE)
require_min_cells_for_pca(seurat_decontx, context = "decontx")
seurat_decontx <- RunPCA(seurat_decontx, seed.use = WORKFLOW_SEED, verbose = FALSE)
seurat_decontx <- RunUMAP(seurat_decontx, dims = 1:30, seed.use = WORKFLOW_SEED)
seurat_decontx <- FindNeighbors(seurat_decontx, dims = 1:30)
seurat_decontx <- FindClusters(seurat_decontx, random.seed = WORKFLOW_SEED)
seurat_decontx <- add_silhouette_to_metadata(seurat_decontx)
saveRDS(seurat_decontx, file = output)
write_cluster_metadata(seurat_decontx, nclusters_output, cluster_ids_output)
