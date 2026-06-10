args <- commandArgs(trailingOnly = TRUE)
filtered <- args[1]
output <- args[2]
nclusters_output <- args[3]
cluster_ids_output <- args[4]

library("Seurat")
library("glmGamPoi")

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

filtered_loaded <- Seurat::Read10X(filtered)
seurat_obj <- CreateSeuratObject(counts = filtered_loaded)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "(?i)^mt-")
seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", seed.use = WORKFLOW_SEED, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, seed.use = WORKFLOW_SEED, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, seed.use = WORKFLOW_SEED)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, random.seed = WORKFLOW_SEED)
seurat_obj <- add_silhouette_to_metadata(seurat_obj)
saveRDS(seurat_obj,file=output)
write_cluster_metadata(seurat_obj, nclusters_output, cluster_ids_output)
