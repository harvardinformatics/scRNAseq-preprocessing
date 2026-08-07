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

seurat <- readRDS(seurat)

seurat_filtered <- subset(seurat, subset = nFeature_RNA > min_nfeature & nCount_RNA > min_ncount & percent.mt < max_mtdna)

# subset() carries the upstream SCT assay, PCA/UMAP embeddings and neighbor graphs into
# the filtered object, but this rule recomputes all of them from RNA counts below. Rebuild
# a minimal counts-only object (preserving percent.mt for the SCTransform regression) so
# that baggage is not held alongside the freshly computed results. Results are unchanged.
seurat_filtered <- CreateSeuratObject(
  counts = GetAssayData(seurat_filtered, assay = "RNA", layer = "counts"),
  meta.data = seurat_filtered@meta.data[, "percent.mt", drop = FALSE]
)
rm(seurat)
gc(verbose = FALSE)

require_min_cells_for_pca(seurat_filtered, context = "posthocfilter_threshold")
seurat_filtered <- SCTransform(seurat_filtered, vars.to.regress = "percent.mt", seed.use = WORKFLOW_SEED, verbose = FALSE)
seurat_filtered <- RunPCA(seurat_filtered, seed.use = WORKFLOW_SEED, verbose = FALSE)
seurat_filtered <- RunUMAP(seurat_filtered, dims = 1:30, seed.use = WORKFLOW_SEED)
seurat_filtered <- FindNeighbors(seurat_filtered, dims = 1:30)
seurat_filtered <- FindClusters(seurat_filtered, random.seed = WORKFLOW_SEED)
seurat_filtered <- add_silhouette_to_metadata(seurat_filtered)
saveRDS(seurat_filtered,file=filtered_output)
write_cluster_metadata(seurat_filtered, nclusters_output, cluster_ids_output)
