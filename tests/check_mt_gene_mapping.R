# Checks that an mt_genes list maps to the row names Seurat gives the listed features.
#
#   check_mt_gene_mapping.R MATRIX_DIR MT_GENES EXPECTED_NAME...
#
# MATRIX_DIR is a 10x matrix directory whose features exercise Seurat's naming (duplicated names that
# Read10X's make.unique suffixes, "_" and "|" that CreateSeuratObject turns into "-"). The object is built as
# tenx2seuratrds.R builds it; the IDs in MT_GENES must map to exactly EXPECTED_NAME..., in row order, and
# percent.mt must be their share of each cell's counts.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("usage: check_mt_gene_mapping.R MATRIX_DIR MT_GENES EXPECTED_NAME...", call. = FALSE)
}
suppressPackageStartupMessages(library(Seurat))
source("workflow/scripts/regression_utils.R")

matrix_dir <- args[1]
expected <- args[-(1:2)]
opts <- parse_script_args(c(paste0("--mt-genes=", args[2]), "--sample=crafted"))
counts <- Read10X(matrix_dir)
obj <- suppressWarnings(CreateSeuratObject(counts = counts))
mapped <- map_mt_genes(opts, read_feature_ids(matrix_dir), obj, context = "check_mt_gene_mapping")
if (!identical(mapped, expected)) {
  stop(sprintf("mapped to %s, expected %s; the object's rows are %s", paste(mapped, collapse = ","),
               paste(expected, collapse = ","), paste(rownames(obj), collapse = ",")), call. = FALSE)
}
obj <- add_percent_mt(obj, mapped, opts, context = "check_mt_gene_mapping")
raw_counts <- LayerData(obj, assay = "RNA", layer = "counts")
share <- 100 * Matrix::colSums(raw_counts[mapped, , drop = FALSE]) / Matrix::colSums(raw_counts)
if (max(abs(obj$percent.mt - share)) > 1e-12) {
  stop("percent.mt is not the listed features' share of the counts", call. = FALSE)
}
if (!identical(obj@misc$percent_mt_features, mapped)) {
  stop("the object does not record the mapped features", call. = FALSE)
}
cat("mapped:", paste(mapped, collapse = ", "), "\n")
