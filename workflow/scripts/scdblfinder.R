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

# scDblFinder fits its doublet classifier in scDblFinder:::.xgbtrain(), once per training round, inside a
# tryCatch() that on error silently keeps the previous score. If every round fails, the calls come from the
# pre-training heuristic score instead of the classifier, and the step would still "succeed". That happens
# with xgboost >= 3.4, which rejects the early_stopping_rounds that scDblFinder < 1.27.12 passes
# (dmlc/xgboost#12370). So watch every training call and stop if one fails on a non-empty training set, or
# if no call returns a model. trace(exit = ) runs as each call returns; it draws no random numbers and
# changes no argument. With SerialParam every call runs in this R process.
# One failure is expected and tolerated: a round that excludes no cells trains on preds[-w, ] with w
# empty, which in R selects zero rows, so xgboost errors and scDblFinder keeps the previous round's trained
# score (seen on the small test-data objects with scDblFinder 1.24.10).
watch_scdblfinder_training <- function() {
  ns <- asNamespace("scDblFinder")
  if (!exists(".xgbtrain", envir = ns, inherits = FALSE)) {
    stop("scdblfinder: scDblFinder ", as.character(packageVersion("scDblFinder")), " has no .xgbtrain(), ",
         "so the classifier-training check in workflow/scripts/scdblfinder.R cannot run; update it.",
         call. = FALSE)
  }
  training_data <- as.name(names(formals(get(".xgbtrain", envir = ns)))[1])
  calls <- new.env(parent = emptyenv())
  calls$trained <- logical(0)
  calls$rows <- integer(0)
  assign(".scdblfinder_training_calls", calls, envir = globalenv())
  invisible(suppressMessages(trace(".xgbtrain", where = ns, print = FALSE, exit = bquote({
    .calls <- get(".scdblfinder_training_calls", envir = globalenv())
    .calls$trained <- c(.calls$trained, !is.null(returnValue(NULL)))
    .calls$rows <- c(.calls$rows, NROW(.(training_data)))
  }))))
  calls
}

stop_unless_scdblfinder_trained <- function(calls) {
  invisible(suppressMessages(untrace(".xgbtrain", where = asNamespace("scDblFinder"))))
  n_calls <- length(calls$trained)
  n_trained <- sum(calls$trained)
  failed_with_data <- sum(!calls$trained & calls$rows > 0L)
  failed_empty <- sum(!calls$trained & calls$rows == 0L)
  if (n_trained == 0L || failed_with_data > 0L) {
    stop(sprintf(paste0(
      "scdblfinder: scDblFinder's classifier training failed (%d of %d training calls returned no model, ",
      "%d of them on a non-empty training set%s). scDblFinder %s silently keeps its previous score when ",
      "training fails, so these would not be the classifier's calls. Known cause: xgboost >= 3.4 ",
      "(installed: %s) with scDblFinder < 1.27.12; keep xgboost < 3.4 in workflow/envs/scrnaseq_r.yml ",
      "and its lockfile."),
      n_calls - n_trained, n_calls, failed_with_data, if (n_calls == 0L) "; no training call was observed" else "",
      as.character(packageVersion("scDblFinder")), as.character(packageVersion("xgboost"))), call. = FALSE)
  }
  if (failed_empty > 0L) {
    message(sprintf(paste0(
      "scdblfinder: classifier trained in %d of %d training calls; the other %d had an empty training set ",
      "(a round that excludes no cells trains on zero rows), so scDblFinder kept the previous round's ",
      "trained score"), n_trained, n_calls, failed_empty))
  } else {
    message(sprintf("scdblfinder: classifier trained in all %d training calls", n_calls))
  }
}

seurat <- readRDS(seurat)
sce <- as.SingleCellExperiment(seurat)
training_calls <- watch_scdblfinder_training()
set.seed(WORKFLOW_SEED)
sce <- scDblFinder(sce, BPPARAM = BiocParallel::SerialParam())
stop_unless_scdblfinder_trained(training_calls)
seurat$scDblFinder.class <- colData(sce)$scDblFinder.class
rm(sce)
gc()
seurat_singlets <- subset(seurat, subset = scDblFinder.class == "singlet")
rm(seurat)
gc()
# percent.mt is kept as metadata only: "(?i)^mt-" finds mitochondrial genes in few species,
# so it is neither filtered on nor regressed out (every species is processed the same way).
seurat_singlets[["percent.mt"]] <- PercentageFeatureSet(seurat_singlets, pattern = "(?i)^mt-")
seurat_singlets <- SCTransform(seurat_singlets, seed.use = WORKFLOW_SEED, verbose = FALSE)
require_min_cells_for_pca(seurat_singlets, context = "scdblfinder")
seurat_singlets <- RunPCA(seurat_singlets, seed.use = WORKFLOW_SEED, verbose = FALSE)
seurat_singlets <- RunUMAP(seurat_singlets, dims = 1:30, seed.use = WORKFLOW_SEED)
seurat_singlets <- FindNeighbors(seurat_singlets, dims = 1:30)
seurat_singlets <- FindClusters(seurat_singlets, random.seed = WORKFLOW_SEED)
seurat_singlets <- add_silhouette_to_metadata(seurat_singlets)
saveRDS(seurat_singlets,file=output)
write_cluster_metadata(seurat_singlets, nclusters_output, cluster_ids_output)
