args <- commandArgs(trailingOnly = TRUE)
seurat <- args[1]
output <- args[2]
markers_output <- args[3]

library("Seurat")
library("DoubletFinder")
library("igraph")
library("tidyverse")
library("tools")
options(future.globals.maxSize = 16 * 1024^3)
seurat <- readRDS(seurat)

sweep.res.list <- paramSweep(seurat, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
optimal_pk <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric),]$pK))

# set value of nExp
nExp_poi <- round(0.15*nrow(seurat@meta.data))
homotypic.prop <- modelHomotypic(seurat$seurat_clusters)
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

seurat <- doubletFinder(seurat, PCs = 1:10, pN = 0.25,
                                 pK = optimal_pk, nExp = nExp_poi.adj,
                                 reuse.pANN = NULL,
                                 sct=TRUE)

classifier_colname <- paste("DF.classifications_0.25",optimal_pk,nExp_poi.adj,sep="_")
seurat_singlets <- subset(seurat, cells = rownames(seurat@meta.data)[seurat@meta.data[[classifier_colname]] == "Singlet"]) 

seurat_singlets[["percent.mt"]] <- PercentageFeatureSet(seurat_singlets, pattern = "(?i)^mt-")
seurat_singlets <- SCTransform(seurat_singlets, vars.to.regress = "percent.mt", verbose = FALSE)
seurat_singlets <- RunPCA(seurat_singlets, verbose = FALSE)
seurat_singlets <- RunUMAP(seurat_singlets, dims = 1:30)
seurat_singlets <- FindNeighbors(seurat_singlets, dims = 1:30)
seurat_singlets <- FindClusters(seurat_singlets)

saveRDS(seurat_singlets,file=output)

markers <- FindAllMarkers(seurat_singlets)
markers$genesymbol <- row.names(markers)
workflow_name <- file_path_sans_ext(basename(markers_output))
workflow_name <- gsub("_markergenes", "", workflow_name)
sig_markers <- as_tibble(markers) %>% filter(p_val_adj<=0.05) %>% mutate(workflow=workflow_name)
write_csv(sig_markers,file=markers_output)
