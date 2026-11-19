# Checks percent.mt and SCTransform's regression in a Seurat object written by a workflow script.
#
#   check_percent_mt.R built   RDS MATRIX_DIR MT_GENES REGRESSED
#       RDS was built from the 10x matrix MATRIX_DIR with the mt_genes list MT_GENES. Its percent.mt and the
#       features it records are recomputed here from the matrix files alone (no Seurat reader).
#   check_percent_mt.R own     RDS FEATURES MT_GENES REGRESSED
#       As built, for objects whose counts are not a 10x matrix's (SoupX, decontX, CellBender): the features
#       come from FEATURES (a 10x matrix directory or a CellBender .h5), percent.mt from the object's counts.
#   check_percent_mt.R carried RDS UPSTREAM_RDS rna|sct|kept REGRESSED
#       RDS was made from UPSTREAM_RDS by a later step, which must reuse the features UPSTREAM_RDS records:
#       "rna" if the step computed percent.mt on the RNA counts, as every step does for a regressed percent.mt;
#       "sct" on the upstream SCT assay (doubletfinder, scdblfinder without regression); "kept" if it kept the
#       upstream values (posthocfilter_threshold without regression). With "rna" and an upstream object whose
#       SCT assay would give other values, the upstream default assay and the gap are reported.
#
# REGRESSED is the comma-separated --vars-to-regress the object was built with, or "none". Its SCTransform
# command must record exactly those variables, and with percent.mt among them every residual (scale.data
# row) must be uncorrelated with percent.mt, which is what the regression does; without, they must not all be.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5 || !args[1] %in% c("built", "own", "carried")) {
  stop("usage: check_percent_mt.R built|own RDS MATRIX_DIR|FEATURES MT_GENES REGRESSED | ",
       "carried RDS UPSTREAM_RDS sct|rna|kept REGRESSED", call. = FALSE)
}
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})
mode <- args[1]
obj <- readRDS(args[2])
regressed <- if (args[5] == "none") NULL else strsplit(args[5], ",", fixed = TRUE)[[1]]

fail <- function(...) stop(sprintf(...), call. = FALSE)
read_gz_lines <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  readLines(con)
}
pct <- obj$percent.mt
if (is.null(pct) || any(!is.finite(pct))) fail("percent.mt is missing or not finite")

if (mode %in% c("built", "own")) {
  source <- args[3]
  lines <- trimws(readLines(args[4]))
  ids <- unique(lines[nzchar(lines) & !startsWith(lines, "#")])
  if (grepl("\\.h5$", source)) {
    feature_ids <- as.character(rhdf5::h5read(source, "matrix/features/id"))
    feature_names <- as.character(rhdf5::h5read(source, "matrix/features/name"))
  } else {
    features <- strsplit(read_gz_lines(file.path(source, "features.tsv.gz")), "\t", fixed = TRUE)
    feature_ids <- vapply(features, `[`, "", 1)
    feature_names <- vapply(features, `[`, "", 2)
  }
  rows <- which(feature_ids %in% ids)
  if (length(rows) == 0) fail("no listed ID is a feature of %s", source)
  # The row names Read10X (or Read_CellBender_h5_Mat) and CreateSeuratObject give these features.
  expected_features <- gsub("[_|]", "-", make.unique(feature_names))[rows]
  if (!identical(obj@misc$percent_mt_features, expected_features)) {
    fail("misc$percent_mt_features is %s, expected %s", paste(obj@misc$percent_mt_features, collapse = ","),
         paste(expected_features, collapse = ","))
  }
  cells <- colnames(obj)
  if (mode == "built") {
    counts <- as(readMM(file.path(source, "matrix.mtx.gz")), "CsparseMatrix")
    colnames(counts) <- read_gz_lines(file.path(source, "barcodes.tsv.gz"))
    expected <- 100 * colSums(counts[rows, cells, drop = FALSE]) / colSums(counts[, cells, drop = FALSE])
  } else {
    counts <- LayerData(obj, assay = "RNA", layer = "counts")
    if (!identical(rownames(counts)[rows], expected_features)) fail("the object's rows are not in feature order")
    expected <- 100 * colSums(counts[rows, cells, drop = FALSE]) / colSums(counts[, cells, drop = FALSE])
  }
} else {
  upstream <- readRDS(args[3])
  stored <- upstream@misc$percent_mt_features
  if (is.null(stored)) fail("the upstream object records no mt_genes features")
  if (!identical(obj@misc$percent_mt_features, stored)) {
    fail("the object does not carry the upstream object's mt_genes features")
  }
  cells <- colnames(obj)
  if (!all(cells %in% colnames(upstream))) fail("cells missing from the upstream object")
  rna_counts <- LayerData(upstream, assay = "RNA", layer = "counts")
  rna_based <- 100 * colSums(rna_counts[stored, cells, drop = FALSE]) / colSums(rna_counts[, cells, drop = FALSE])
  sct_based <- NULL
  if ("SCT" %in% Assays(upstream)) {
    sct_counts <- LayerData(upstream, assay = "SCT", layer = "counts")
    used <- intersect(stored, rownames(sct_counts))
    sct_based <- 100 * colSums(sct_counts[used, cells, drop = FALSE]) / upstream$nCount_SCT[cells]
  }
  expected <- switch(args[4], rna = rna_based, sct = sct_based, kept = upstream$percent.mt[cells])
  if (is.null(expected)) fail("the upstream object has no SCT assay")
  if (args[4] == "rna" && !is.null(sct_based)) {
    # The check can tell RNA from SCT only if they differ.
    gap <- max(abs(rna_based - sct_based))
    if (gap < 1e-6) fail("RNA- and SCT-based percent.mt coincide, so the RNA assay's use cannot be checked")
    cat(sprintf("upstream default assay: %s; RNA- and SCT-based percent.mt differ by up to %.3g\n",
                DefaultAssay(upstream), gap))
  }
}
diff <- max(abs(pct[cells] - expected))
if (diff > 1e-9) fail("percent.mt differs from the expected values by up to %g", diff)
if (max(pct) <= 0) fail("percent.mt is 0 in every cell")

sct_command <- obj@commands[["SCTransform.RNA"]]
if (is.null(sct_command)) fail("no SCTransform.RNA command recorded")
if (!identical(sct_command@params$vars.to.regress, regressed)) {
  fail("SCTransform's vars.to.regress was %s, expected %s", deparse(sct_command@params$vars.to.regress),
       deparse(regressed))
}
scale_data <- LayerData(obj, assay = "SCT", layer = "scale.data")
varying <- apply(scale_data, 1, sd) > 0
cors <- abs(apply(scale_data[varying, , drop = FALSE], 1, function(row) cor(row, pct[colnames(scale_data)])))
if ("percent.mt" %in% regressed) {
  if (max(cors) > 1e-8) fail("percent.mt was regressed out, but a residual correlates with it (|r| up to %g)", max(cors))
} else if (median(cors) < 1e-3) {
  fail("percent.mt was not regressed out, yet no residual correlates with it")
}
cat(sprintf("%s: percent.mt matches (max diff %.2g) from %d features; vars.to.regress = %s; max |cor(residual, percent.mt)| = %.2g\n",
            args[2], diff, length(obj@misc$percent_mt_features), deparse(regressed), max(cors)))
