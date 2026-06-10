args <- commandArgs(trailingOnly = TRUE)
seurat <- args[1]
output <- args[2]
nclusters_output <- args[3]
cluster_ids_output <- args[4]

library("Seurat")
library("scDblFinder")

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
sce <- as.SingleCellExperiment(seurat)
set.seed(WORKFLOW_SEED)
sce <- scDblFinder(sce, BPPARAM = BiocParallel::SerialParam())
seurat$scDblFinder.class <- colData(sce)$scDblFinder.class
rm(sce)
gc()
seurat_singlets <- subset(seurat, subset = scDblFinder.class == "singlet")
rm(seurat)
gc()
seurat_singlets[["percent.mt"]] <- PercentageFeatureSet(seurat_singlets, pattern = "(?i)^mt-")
seurat_singlets <- SCTransform(seurat_singlets, vars.to.regress = "percent.mt", seed.use = WORKFLOW_SEED, verbose = FALSE)
seurat_singlets <- RunPCA(seurat_singlets, seed.use = WORKFLOW_SEED, verbose = FALSE)
seurat_singlets <- RunUMAP(seurat_singlets, dims = 1:30, seed.use = WORKFLOW_SEED)
seurat_singlets <- FindNeighbors(seurat_singlets, dims = 1:30)
seurat_singlets <- FindClusters(seurat_singlets, random.seed = WORKFLOW_SEED)
seurat_singlets <- add_silhouette_to_metadata(seurat_singlets)
saveRDS(seurat_singlets,file=output)
write_cluster_metadata(seurat_singlets, nclusters_output, cluster_ids_output)
