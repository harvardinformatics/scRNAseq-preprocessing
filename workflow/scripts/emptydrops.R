args <- commandArgs(trailingOnly = TRUE)
raw <- args[1]
output <- args[2]

library("Seurat")
library("DropletUtils")
library("scater")

droputil_rawdata <-read10xCounts(raw)
colnames(droputil_rawdata) <- droputil_rawdata$Barcode
rownames(droputil_rawdata) <- make.unique(rowData(droputil_rawdata)$Symbol)
dropletutils_out <-emptyDrops(counts(droputil_rawdata))
is.cell <-!is.na(dropletutils_out$FDR) & dropletutils_out$FDR<0.01
cell_barcodes <- colnames(droputil_rawdata)[which(is.cell)]
droputil_filtered <- droputil_rawdata[,is.cell]
droputil_filtered <- logNormCounts(droputil_filtered)
seurat_droputil_filtered <- as.Seurat(droputil_filtered, counts = "counts",data="logcounts") 
seurat_droputil_filtered[["percent.mt"]] <- PercentageFeatureSet(seurat_droputil_filtered, pattern = "(?i)^mt-")
saveRDS(seurat_droputil_filtered,file=output)
