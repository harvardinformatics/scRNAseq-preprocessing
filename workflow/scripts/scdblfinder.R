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

# scDblFinder fits its doublet classifier in scDblFinder:::.xgbtrain(), inside a tryCatch() that on
# error silently keeps the pre-training score, so the step would "succeed" with heuristic calls
# instead of the classifier's. That happens with xgboost >= 3.4, which rejects the
# early_stopping_rounds that scDblFinder < 1.27.12 passes (dmlc/xgboost#12370). Count the training
# calls that return a model and stop if any fails. trace(exit = ) runs as each call returns; it draws
# no random numbers and changes no argument. With SerialParam every call runs in this R process.
watch_scdblfinder_training <- function() {
  ns <- asNamespace("scDblFinder")
  if (!exists(".xgbtrain", envir = ns, inherits = FALSE)) {
    stop("scdblfinder: scDblFinder ", as.character(packageVersion("scDblFinder")), " has no .xgbtrain(), ",
         "so the classifier-training check in workflow/scripts/scdblfinder.R cannot run; update it.",
         call. = FALSE)
  }
  calls <- new.env(parent = emptyenv())
  calls$trained <- 0L
  calls$failed <- 0L
  assign(".scdblfinder_training_calls", calls, envir = globalenv())
  invisible(suppressMessages(trace(".xgbtrain", where = ns, print = FALSE, exit = quote({
    .calls <- get(".scdblfinder_training_calls", envir = globalenv())
    .outcome <- if (is.null(returnValue(NULL))) "failed" else "trained"
    assign(.outcome, get(.outcome, envir = .calls) + 1L, envir = .calls)
  }))))
  calls
}

stop_unless_scdblfinder_trained <- function(calls) {
  invisible(suppressMessages(untrace(".xgbtrain", where = asNamespace("scDblFinder"))))
  n_calls <- calls$trained + calls$failed
  if (n_calls == 0L || calls$failed > 0L) {
    stop(sprintf(paste0(
      "scdblfinder: scDblFinder's doublet classifier did not train (%d of %d training calls failed%s). ",
      "scDblFinder %s silently falls back to its pre-training score, so these would not be scDblFinder's ",
      "calls. Known cause: xgboost >= 3.4 (installed: %s) with scDblFinder < 1.27.12; keep xgboost < 3.4 ",
      "in workflow/envs/scrnaseq_r.yml and its lockfile."),
      calls$failed, n_calls, if (n_calls == 0L) "; no training call was observed" else "",
      as.character(packageVersion("scDblFinder")), as.character(packageVersion("xgboost"))), call. = FALSE)
  }
  message(sprintf("scdblfinder: classifier trained in all %d training calls", n_calls))
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
