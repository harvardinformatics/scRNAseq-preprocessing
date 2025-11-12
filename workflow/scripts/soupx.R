args <- commandArgs(trailingOnly = TRUE)
filtered <- args[1]
raw <- args[2]
seurat_base <- args[3]
output <- args[4]

library("Seurat")
library("glmGamPoi")
library("SoupX")

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
seurat_soupx[["percent.mt"]] <- PercentageFeatureSet(seurat_soupx, pattern = "(?i)^mt-")
seurat_soupx <- SCTransform(seurat_soupx, vars.to.regress = "percent.mt", verbose = FALSE)
seurat_soupx <- RunPCA(seurat_soupx, verbose = FALSE)
seurat_soupx <- RunUMAP(seurat_soupx, dims = 1:30)
seurat_soupx <- FindNeighbors(seurat_soupx, dims = 1:30)
seurat_soupx <- FindClusters(seurat_soupx)
saveRDS(seurat_soupx,file=output)
