# SCALABILITY NOTE: doubletFinder() augments the full dataset with ~25% synthetic doublets
# and then builds a DENSE pairwise distance matrix (fields::rdist) over all cells. That matrix
# is O(N^2): ~78 GB at 77k cells, ~132 GB at 100k, ~298 GB at 150k, ~530 GB at 200k. paramSweep
# is not the driver (it caps at a 10k-cell subsample); the main classification step is. This is
# inherent to the DoubletFinder algorithm and cannot be tuned away here, so for large datasets
# (roughly >100k cells, e.g. emptyDrops outputs) prefer scDblFinder (the workflow's other
# doublet_removal_method), which does not materialize a full distance matrix and scales far
# better. Memory for this rule is provisioned by input size in workflow/rules/doubletfinder.smk.

args <- commandArgs(trailingOnly = TRUE)
seurat <- args[1]
output <- args[2]
nclusters_output <- args[3]
cluster_ids_output <- args[4]

library("Seurat")
library("DoubletFinder")
library("igraph")
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

set.seed(WORKFLOW_SEED)
sweep.res.list <- paramSweep(seurat, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
optimal_pk <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric),]$pK))
rm(sweep.res.list, sweep.stats, bcmvn)
gc()

# set value of nExp
nExp_poi <- round(0.15*nrow(seurat@meta.data))
homotypic.prop <- modelHomotypic(seurat$seurat_clusters)
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

set.seed(WORKFLOW_SEED)
seurat <- doubletFinder(seurat, PCs = 1:10, pN = 0.25,
                                 pK = optimal_pk, nExp = nExp_poi.adj,
                                 reuse.pANN = NULL,
                                 sct=TRUE)

classifier_colname <- paste("DF.classifications_0.25",optimal_pk,nExp_poi.adj,sep="_")
seurat_singlets <- subset(seurat, cells = rownames(seurat@meta.data)[seurat@meta.data[[classifier_colname]] == "Singlet"]) 
rm(seurat)
gc()

seurat_singlets[["percent.mt"]] <- PercentageFeatureSet(seurat_singlets, pattern = "(?i)^mt-")
seurat_singlets <- SCTransform(seurat_singlets, vars.to.regress = "percent.mt", seed.use = WORKFLOW_SEED, verbose = FALSE)
require_min_cells_for_pca(seurat_singlets, context = "doubletfinder")
seurat_singlets <- RunPCA(seurat_singlets, seed.use = WORKFLOW_SEED, verbose = FALSE)
seurat_singlets <- RunUMAP(seurat_singlets, dims = 1:30, seed.use = WORKFLOW_SEED)
seurat_singlets <- FindNeighbors(seurat_singlets, dims = 1:30)
seurat_singlets <- FindClusters(seurat_singlets, random.seed = WORKFLOW_SEED)
seurat_singlets <- add_silhouette_to_metadata(seurat_singlets)

saveRDS(seurat_singlets,file=output)
write_cluster_metadata(seurat_singlets, nclusters_output, cluster_ids_output)
