args <- commandArgs(trailingOnly = TRUE)
cellbender_h5 <- args[1]
output <- args[2]
nclusters_output <- args[3]
cluster_ids_output <- args[4]

library("Seurat")
library("scCustomize")
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

mat <- Read_CellBender_h5_Mat(cellbender_h5)
# CellBender's *_filtered.h5 can still contain a rare barcode whose counts were
# all subtracted to zero during background removal. A zero-count cell yields
# log_umi = log10(0) = -Inf, which makes SCTransform's make_cell_attr abort.
# min.features = 1 drops only these empty barcodes and leaves all real cells intact.
seurat <- CreateSeuratObject(mat, min.features = 1)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "(?i)^mt-")
seurat <- SCTransform(seurat, vars.to.regress = "percent.mt", seed.use = WORKFLOW_SEED, verbose = FALSE)
require_min_cells_for_pca(seurat, context = "cellbender2seurat")
seurat <- RunPCA(seurat, seed.use = WORKFLOW_SEED, verbose = FALSE)
seurat <- RunUMAP(seurat, dims = 1:30, seed.use = WORKFLOW_SEED)
seurat <- FindNeighbors(seurat, dims = 1:30)
seurat <- FindClusters(seurat, random.seed = WORKFLOW_SEED)
seurat <- add_silhouette_to_metadata(seurat)
saveRDS(seurat,file=output)
write_cluster_metadata(seurat, nclusters_output, cluster_ids_output)
