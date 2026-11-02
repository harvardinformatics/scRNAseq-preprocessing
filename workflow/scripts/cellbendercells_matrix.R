# Raw counts of the barcodes CellBender calls as cells, without its background removal.
# CellBender calls cells and removes ambient RNA in one model, so the cellbender_fromraw object
# carries both effects; this matrix carries only the calls. It is written as a 10x matrix
# directory and built by tenx2seuratrds.R (rule cellbendercells2seurat) exactly like the
# CellRanger filtered matrix, so it compares directly with the CellRanger and emptyDrops
# cell-calling objects.
args <- commandArgs(trailingOnly = TRUE)
raw <- args[1]
cellbender_h5 <- args[2]
matrix_outdir <- args[3]

library("DropletUtils")
library("rhdf5")

# The *_filtered.h5 holds exactly the barcodes CellBender called as cells; its counts are
# background-removed and are not used here.
cell_barcodes <- as.character(h5read(cellbender_h5, "matrix/barcodes"))
if (length(cell_barcodes) == 0) {
    stop("no cell barcodes in ", cellbender_h5)
}
if (anyDuplicated(cell_barcodes)) {
    stop("duplicate cell barcodes in ", cellbender_h5)
}

# Same reader as emptydrops.R: the counts are the raw counts every cell caller starts from.
raw_counts <- read10xCounts(raw)
raw_barcodes <- raw_counts$Barcode
missing <- setdiff(cell_barcodes, raw_barcodes)
if (length(missing) > 0) {
    stop(length(missing), " of ", length(cell_barcodes), " CellBender cell barcodes are not in the ",
         "raw matrix (e.g. ", paste(head(missing, 3), collapse = ", "), "); the raw matrix and the ",
         "raw h5 CellBender read do not match")
}

# Keep the raw matrix's barcode order, as CellRanger's filtered matrix does.
keep <- which(raw_barcodes %in% cell_barcodes)
cm <- counts(raw_counts)[, keep, drop = FALSE]
if (!inherits(cm, "dgCMatrix")) {
    cm <- as(cm, "dgCMatrix")
}
barcodes <- raw_barcodes[keep]
umis <- Matrix::colSums(cm)
cat(sprintf("raw barcodes: %d\nCellBender cell calls: %d\nraw UMIs per called cell: min %g, median %g, max %g\n",
            length(raw_barcodes), length(keep), min(umis), median(umis), max(umis)))

# A barcode with no counts makes SCTransform abort (log10(0)), as in cellbender2seurat.R. CellBender
# calls cells from these raw counts, so none is expected; drop any rather than fail the build.
empty <- umis == 0
if (any(empty)) {
    cat(sprintf("dropping %d called barcode(s) with no raw counts\n", sum(empty)))
    cm <- cm[, !empty, drop = FALSE]
    barcodes <- barcodes[!empty]
}

features <- rowData(raw_counts)
gene_type <- if ("Type" %in% colnames(features)) features$Type else "Gene Expression"
write10xCounts(matrix_outdir, cm, barcodes = barcodes, gene.id = features$ID,
               gene.symbol = features$Symbol, gene.type = gene_type, version = "3", overwrite = TRUE)
cat(sprintf("wrote %d cells x %d features to %s\n", ncol(cm), nrow(cm), matrix_outdir))
