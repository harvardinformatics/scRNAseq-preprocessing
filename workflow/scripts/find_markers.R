args <- commandArgs(trailingOnly = TRUE)
rds_input <- args[1]
cluster_id <- args[2]
output_path <- args[3]

library("Seurat")
library("tidyverse")
library("tools")

options(future.globals.maxSize = 16 * 1024^3)

seurat_obj <- readRDS(rds_input)
Idents(seurat_obj) <- "seurat_clusters"

markers <- FindMarkers(seurat_obj, ident.1 = cluster_id)
markers_tbl <- markers %>%
  as.data.frame() %>%
  rownames_to_column(var = "genesymbol") %>%
  as_tibble()

workflow_name <- file_path_sans_ext(basename(output_path))
workflow_name <- sub("_markergenes_cluster.*$", "", workflow_name)

sig_markers <- markers_tbl %>%
  filter(p_val_adj <= 0.05) %>%
  mutate(cluster = cluster_id, workflow = workflow_name)

write_csv(sig_markers, file = output_path)
