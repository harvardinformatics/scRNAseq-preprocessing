args <- commandArgs(trailingOnly = TRUE)
rds_input <- args[1]
cluster_id <- args[2]
output_path <- args[3]

library("Seurat")
library("tidyverse")
library("tools")

options(future.globals.maxSize = 16 * 1024^3)

source("workflow/scripts/silhouette_utils.R")
WORKFLOW_SEED <- set_workflow_seed()

seurat_obj <- readRDS(rds_input)
Idents(seurat_obj) <- "seurat_clusters"

workflow_name <- file_path_sans_ext(basename(output_path))
workflow_name <- sub("_markergenes_cluster.*$", "", workflow_name)

# FindMarkers (via ValidateCellGroups) requires >= 3 cells in the cluster.
# Clustering can legitimately emit micro-clusters (outliers / residual doublets)
# with 1-2 cells; those cannot yield markers. Skip them gracefully by writing a
# schema-correct empty file so combine_markers can still concatenate all clusters.
n_cells_in_cluster <- sum(Idents(seurat_obj) == cluster_id)
empty_markers <- tibble(
  genesymbol = character(),
  p_val = numeric(),
  avg_log2FC = numeric(),
  pct.1 = numeric(),
  pct.2 = numeric(),
  p_val_adj = numeric(),
  cluster = character(),
  workflow = character()
)

if (n_cells_in_cluster < 3) {
  message(sprintf(
    "Cluster %s has %d cell(s) (< 3); skipping marker detection and writing empty output.",
    cluster_id, n_cells_in_cluster
  ))
  write_csv(empty_markers, file = output_path)
} else {
  set.seed(WORKFLOW_SEED)
  markers <- FindMarkers(seurat_obj, ident.1 = cluster_id)
  markers_tbl <- markers %>%
    as.data.frame() %>%
    rownames_to_column(var = "genesymbol") %>%
    as_tibble()

  sig_markers <- markers_tbl %>%
    filter(p_val_adj <= 0.05) %>%
    mutate(cluster = cluster_id, workflow = workflow_name)

  write_csv(sig_markers, file = output_path)
}
