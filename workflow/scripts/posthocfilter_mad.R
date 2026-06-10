args <- commandArgs(trailingOnly = TRUE)
seurat <- args[1]
filtered_output <- args[2]
nclusters_output <- args[3]
cluster_ids_output <- args[4]

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
sce <-as.SingleCellExperiment(seurat)
mito_genes <- grep("(?i)^mt-", rownames(sce), value = TRUE)
qc <- perCellQCMetrics(sce, subsets = list(Mito = mito_genes))
high_mito <- isOutlier(qc$subsets_Mito_percent, nmads=3, type="higher")
low_umi     <- isOutlier(qc$sum, nmads = 3, type = "lower")
low_feature <- isOutlier(qc$detected, nmads = 3, type = "lower")
discard <- high_mito | low_umi | low_feature
cells_to_keep <- colnames(seurat)[!discard]
seurat_filtered<- subset(seurat, cells = cells_to_keep)

seurat_filtered[["percent.mt"]] <- PercentageFeatureSet(seurat_filtered, pattern = "(?i)^mt-")
seurat_filtered <- SCTransform(seurat_filtered, vars.to.regress = "percent.mt", seed.use = WORKFLOW_SEED, verbose = FALSE)
seurat_filtered <- RunPCA(seurat_filtered, seed.use = WORKFLOW_SEED, verbose = FALSE)
seurat_filtered <- RunUMAP(seurat_filtered, dims = 1:30, seed.use = WORKFLOW_SEED)
seurat_filtered <- FindNeighbors(seurat_filtered, dims = 1:30)
seurat_filtered <- FindClusters(seurat_filtered, random.seed = WORKFLOW_SEED)
seurat_filtered <- add_silhouette_to_metadata(seurat_filtered)
saveRDS(seurat_filtered,file=filtered_output)
write_cluster_metadata(seurat_filtered, nclusters_output, cluster_ids_output)
