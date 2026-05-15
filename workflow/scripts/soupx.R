args <- commandArgs(trailingOnly = TRUE)
filtered <- args[1]
raw <- args[2]
seurat_base <- args[3]
output <- args[4]
nclusters_output <- args[5]
cluster_ids_output <- args[6]

library("Seurat")
library("glmGamPoi")
library("SoupX")

options(future.globals.maxSize = 16 * 1024^3)

source("workflow/scripts/silhouette_utils.R")

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
raw_matrix <- Seurat::Read10X(raw)


seurat_base <- readRDS(seurat_base)
soup_channel <- SoupX::SoupChannel(tod = raw_matrix,toc=filtered_matrix,
                is10X = TRUE)
soup_channel$tod <- raw
soup_channel <- SoupX::setClusters(soup_channel, 
                                   clusters = as.factor(Idents(seurat_base)))
soup_channel <- setDR(soup_channel, 
                DR=Seurat::Embeddings(seurat_base, "umap"))
soup_channel <- autoEstCont(soup_channel)
corrected_counts <- adjustCounts(soup_channel,roundToInt=TRUE)
seurat_soupx <- CreateSeuratObject(counts = corrected_counts)
rm(filtered_matrix, raw_matrix, seurat_base, soup_channel, corrected_counts)
gc()
seurat_soupx[["percent.mt"]] <- PercentageFeatureSet(seurat_soupx, pattern = "(?i)^mt-")
seurat_soupx <- SCTransform(seurat_soupx, vars.to.regress = "percent.mt", verbose = FALSE)
seurat_soupx <- RunPCA(seurat_soupx, verbose = FALSE)
seurat_soupx <- RunUMAP(seurat_soupx, dims = 1:30)
seurat_soupx <- FindNeighbors(seurat_soupx, dims = 1:30)
seurat_soupx <- FindClusters(seurat_soupx)
seurat_soupx <- add_silhouette_to_metadata(seurat_soupx)
saveRDS(seurat_soupx,file=output)
write_cluster_metadata(seurat_soupx, nclusters_output, cluster_ids_output)
