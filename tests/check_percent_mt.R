# Checks percent.mt and SCTransform's regression in a Seurat object written by a workflow script. In every saved
# object, percent.mt must be the share of each cell's RNA counts in its mitochondrial features: those of the
# sample's mt_genes list or, without one, those whose names match "(?i)^mt-".
#
#   check_percent_mt.R built   RDS MATRIX_DIR MT_GENES REGRESSED
#       RDS was built from the 10x matrix MATRIX_DIR with the mt_genes list MT_GENES, or with the pattern if
#       MT_GENES is "pattern". Its percent.mt and the features it records are recomputed here from the matrix
#       files alone (no Seurat reader).
#   check_percent_mt.R own     RDS FEATURES MT_GENES REGRESSED
#       As built, for objects whose counts are not a 10x matrix's (SoupX, decontX, CellBender): the features
#       come from FEATURES (a 10x matrix directory or a CellBender .h5), percent.mt from the object's RNA counts.
#   check_percent_mt.R carried RDS UPSTREAM_RDS REGRESSED
#       RDS was made from UPSTREAM_RDS by a later step, which must reuse the features UPSTREAM_RDS records (the
#       pattern's if it records none) and compute percent.mt on the RNA counts, whatever the upstream object's
#       default assay and whatever percent.mt it carries. If UPSTREAM_RDS has an SCT assay, the percent.mt
#       computed on it must differ, or the check could not tell the assays apart; the upstream default assay
#       and the gap are reported.
#
# REGRESSED is the comma-separated --vars-to-regress the object was built with, or "none". Its SCTransform
# command must record exactly those variables, and with percent.mt among them every residual (scale.data
# row) must be uncorrelated with percent.mt, which is what the regression does; without, they must not all be.
MT_PATTERN <- "(?i)^mt-"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || !args[1] %in% c("built", "own", "carried") ||
    length(args) != if (args[1] == "carried") 4 else 5) {
  stop("usage: check_percent_mt.R built|own RDS MATRIX_DIR|FEATURES MT_GENES|pattern REGRESSED | ",
       "carried RDS UPSTREAM_RDS REGRESSED", call. = FALSE)
}
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})
mode <- args[1]
obj <- readRDS(args[2])
regressed <- if (args[length(args)] == "none") NULL else strsplit(args[length(args)], ",", fixed = TRUE)[[1]]

fail <- function(...) stop(sprintf(...), call. = FALSE)
read_gz_lines <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  readLines(con)
}
matching_pattern <- function(names) grep(MT_PATTERN, names, perl = TRUE)
rna_share <- function(counts, rows, cells) {
  100 * colSums(counts[rows, cells, drop = FALSE]) / colSums(counts[, cells, drop = FALSE])
}
pct <- obj$percent.mt
if (is.null(pct) || any(!is.finite(pct))) fail("percent.mt is missing or not finite")
cells <- colnames(obj)

if (mode %in% c("built", "own")) {
  source <- args[3]
  if (grepl("\\.h5$", source)) {
    feature_ids <- as.character(rhdf5::h5read(source, "matrix/features/id"))
    feature_names <- as.character(rhdf5::h5read(source, "matrix/features/name"))
  } else {
    features <- strsplit(read_gz_lines(file.path(source, "features.tsv.gz")), "\t", fixed = TRUE)
    feature_ids <- vapply(features, `[`, "", 1)
    feature_names <- vapply(features, `[`, "", 2)
  }
  # The row names Read10X (or Read_CellBender_h5_Mat) and CreateSeuratObject give these features.
  row_names <- gsub("[_|]", "-", make.unique(feature_names))
  if (args[4] == "pattern") {
    rows <- matching_pattern(row_names)
    if (length(rows) == 0) fail("no feature name of %s matches %s", source, MT_PATTERN)
    recorded <- NULL
  } else {
    lines <- trimws(readLines(args[4]))
    ids <- unique(lines[nzchar(lines) & !startsWith(lines, "#")])
    rows <- which(feature_ids %in% ids)
    if (length(rows) == 0) fail("no listed ID is a feature of %s", source)
    recorded <- row_names[rows]
  }
  if (!identical(obj@misc$percent_mt_features, recorded)) {
    fail("misc$percent_mt_features is %s, expected %s", deparse(obj@misc$percent_mt_features), deparse(recorded))
  }
  if (mode == "built") {
    counts <- as(readMM(file.path(source, "matrix.mtx.gz")), "CsparseMatrix")
    colnames(counts) <- read_gz_lines(file.path(source, "barcodes.tsv.gz"))
  } else {
    counts <- LayerData(obj, assay = "RNA", layer = "counts")
    if (!identical(rownames(counts), row_names)) fail("the object's rows are not the features of %s", source)
  }
  expected <- rna_share(counts, rows, cells)
  n_features <- length(rows)
} else {
  upstream <- readRDS(args[3])
  recorded <- upstream@misc$percent_mt_features
  if (!identical(obj@misc$percent_mt_features, recorded)) {
    fail("the object does not carry the upstream object's mt_genes features")
  }
  if (!all(cells %in% colnames(upstream))) fail("cells missing from the upstream object")
  mt_rows <- function(counts) {
    if (is.null(recorded)) matching_pattern(rownames(counts)) else which(rownames(counts) %in% recorded)
  }
  rna_counts <- LayerData(upstream, assay = "RNA", layer = "counts")
  rows <- mt_rows(rna_counts)
  if (length(rows) == 0) fail("the upstream object's RNA assay has no mitochondrial feature")
  expected <- rna_share(rna_counts, rows, cells)
  n_features <- length(rows)
  if ("SCT" %in% Assays(upstream)) {
    sct_counts <- LayerData(upstream, assay = "SCT", layer = "counts")
    gap <- max(abs(expected - rna_share(sct_counts, mt_rows(sct_counts), cells)))
    if (gap < 1e-6) fail("RNA- and SCT-based percent.mt coincide, so the RNA assay's use cannot be checked")
    cat(sprintf("upstream default assay: %s; RNA- and SCT-based percent.mt differ by up to %.3g\n",
                DefaultAssay(upstream), gap))
  }
}
diff <- max(abs(pct[cells] - expected))
if (diff > 1e-9) fail("percent.mt differs from the RNA-based values by up to %g", diff)
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
cat(sprintf("%s: percent.mt matches the RNA-based values (max diff %.2g) from %d %s features; vars.to.regress = %s; max |cor(residual, percent.mt)| = %.2g\n",
            args[2], diff, n_features, if (is.null(recorded)) "pattern-matched" else "listed", deparse(regressed),
            max(cors)))
