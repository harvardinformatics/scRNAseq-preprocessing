library(tidyverse)
library(Seurat)
library(Matrix)
library(cluster)
library(purrr)

### File import/export functions ###

#' @title
#' Export Seurat object to matrix and associated files
#' 
#' @description
#' export_seurat_to_files takes a Seurat scRNA-seq object, and exports four files:
#'  * a Matrix Market format expression matrix
#'  * a feature table, which for a standard 10x library consists of gene symbols
#'  * a table of cell barcodes
#'  * a tabular extraction of the Seurat object metadata
#'  These four files are all required inputs for cell type annotation with CellTypist.
#'  
#' @param obj        Seurat object
#' @param prefix     Character prefix for output files (no extension)
#' @param assay      Assay to use (default "RNA")
#' @param layer      Layer to use for v5 objects (default "counts"; set to NULL for v3-style)
#' @return Invisibly returns \code{NULL}. Files are written to disk as a side effect.
export_seurat_to_files <- function(obj,
                                   prefix,
                                   assay = "RNA",
                                   layer = "counts") {
  # Set default assay
  DefaultAssay(obj) <- assay
  
  # Get expression matrix
  if (is.null(layer)) {
    mat <- GetAssayData(obj, assay = assay, slot = "counts")
  } else {
    mat <- GetAssayData(obj, assay = assay, layer = layer)
  }
  
  # Coerce to standard dgCMatrix
  mat <- as(mat, "dgCMatrix")
  
  # Build filenames
  mtx_file       <- paste0(prefix, ".mtx")
  genes_file     <- paste0("genes_", prefix, ".tsv")
  barcodes_file  <- paste0("barcodes_", prefix, ".tsv")
  metadata_file  <- paste0("metadata_", prefix, ".csv")
  
  # Write matrix
  Matrix::writeMM(mat, file = mtx_file)
  
  # Write genes
  write.table(
    data.frame(gene = rownames(mat)),
    file      = genes_file,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  # Write barcodes
  write.table(
    data.frame(barcode = colnames(mat)),
    file      = barcodes_file,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  # Write metadata
  meta <- obj@meta.data
  write.csv(
    meta,
    file      = metadata_file,
    row.names = TRUE
  )
}


### Adding metadata to Seurat objects ###

#' @title
#' Add silhouette widths to Seurat object metadata
#'
#' @description
#' Computes silhouette widths based on a chosen dimensional reduction (e.g. PCA)
#' and a clustering column, then stores the per-cell silhouette values in
#' \code{seurat_obj@meta.data}.
#'
#' @param seurat_obj A \code{Seurat} object containing cell embeddings and metadata.
#' @param cluster_col Character string; name of the metadata column containing
#'   cluster assignments (default: \code{"seurat_clusters"}).
#' @param reduction Character string; name of the reduction to use for distance
#'   calculation (e.g. \code{"pca"}, \code{"umap"}). Default is \code{"pca"}.
#' @param dims Integer vector; dimensions of the chosen reduction to use when
#'   computing distances (default: \code{1:30}).
#' @param sil_col_name Character string; name of the new metadata column in which
#'   silhouette widths will be stored (default: \code{"silhouette"}).
#'
#' @return A \code{Seurat} object with an additional column in \code{meta.data}
#'   containing per-cell silhouette widths.
#'
#' @details
#' Cluster labels in \code{cluster_col} are converted to integer factors before
#' calling \code{cluster::silhouette}. Distances are computed using
#' \code{stats::dist} on the selected embedding dimensions.
#'
#' @importFrom Seurat Embeddings
#' @importFrom cluster silhouette
#' @importFrom stats dist
add_silhouette_to_metadata <- function(
    seurat_obj,
    cluster_col   = "seurat_clusters",
    reduction     = "pca",
    dims          = 1:30,
    sil_col_name  = "silhouette_width"
) {
  # Get cluster labels as integer vector
  if (! cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", cluster_col, "not found in meta.data"))
  }
  clust <- seurat_obj@meta.data[[cluster_col]]
  # silhouette() expects integers
  clust_int <- as.integer(as.factor(clust))
  
  # Get coordinates for distance calculation (e.g., PCA embeddings)
  emb <- Embeddings(seurat_obj, reduction = reduction)[, dims, drop = FALSE]
  
  # Compute Euclidean distance matrix
  d <- dist(emb)   # or use another distance if you prefer
  
  # Compute silhouette
  sil <- silhouette(clust_int, d)
  
  # sil is in the same cell order as emb/meta.data
  sil_values <- sil[, "sil_width"]
  
  # Add to meta.data
  seurat_obj@meta.data[[sil_col_name]] <- sil_values
  
  return(seurat_obj)
}


### Functions for analyzing cell barcode composition

#' @title
#' Compute Jaccard similarity between two sets
#'
#' @description
#' Calculates the Jaccard similarity coefficient between two vectors treated
#' as sets, defined as the size of their intersection divided by the size of
#' their union.
#'
#' @param set1 A vector representing the first set (will be treated as unique elements).
#' @param set2 A vector representing the second set (will be treated as unique elements).
#'
#' @return A numeric value between 0 and 1 giving the Jaccard similarity:
#'   \deqn{|set1 ∩ set2| / |set1 ∪ set2|}.
#'
#' @examples
#' jaccard_similarity(c("a", "b", "c"), c("b", "c", "d"))
#' jaccard_similarity(1:3, 2:4)
jaccard_similarity <- function(set1, set2) {
  intersect_length <- length(intersect(set1, set2))
  union_length <- length(set1) + length(set2) - intersect_length
  intersect_length / union_length
}

#' @title
#' Plot a Jaccard similarity heatmap between two clusterings
#'
#' @description
#' Computes pairwise Jaccard similarity between clusters from two different
#' metadata tables (typically corresponding to two Seurat objects / clusterings)
#' and plots the result as a heatmap.
#'
#' @param meta1 A data frame (e.g. \code{seurat_obj@meta.data}) containing
#'   metadata for the first clustering. Row names must be cell barcodes and
#'   there must be a \code{seurat_clusters} column.
#' @param meta2 A data frame containing metadata for the second clustering.
#'   Row names must be cell barcodes and there must be a \code{seurat_clusters}
#'   column.
#' @param name1 Character string; label to use for the x-axis (e.g. name of the
#'   first clustering method).
#' @param name2 Character string; label to use for the y-axis (e.g. name of the
#'   second clustering method).
#' @param threshold Numeric value in \code{[0, 1]} used to highlight cells in the
#'   heatmap: cluster pairs with Jaccard similarity greater than or equal to this
#'   threshold are marked with a star symbol.
#'
#' @return A \code{ggplot2} object representing the Jaccard similarity heatmap.
#'
#' @details
#' Cell barcodes are taken from the row names of \code{meta1} and \code{meta2}.
#' The function performs a full join on barcodes, then splits rows by cluster
#' identity from each metadata table and computes the Jaccard similarity between
#' all pairs of clusters using \code{\link{jaccard_similarity}}.
#'
#' @importFrom tibble tibble
#' @importFrom dplyr full_join
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient geom_point
#'   scale_x_discrete scale_y_discrete labs theme_minimal theme element_blank
#'   element_text
#'
#' @examples
#' \dontrun{
#' p <- jaccard_heatmap_plot(
#'   meta1      = seurat_obj1@meta.data,
#'   meta2      = seurat_obj2@meta.data,
#'   name1      = "Clustering 1",
#'   name2      = "Clustering 2",
#'   threshold  = 0.4
#' )
#' print(p)
#' }
jaccard_heatmap_plot <- function(meta1, meta2,
                                 clusterid1, clusterid2,
                                 name1, name2, threshold) {
  
  clusters1 <- tibble::tibble(
    cellbarcode = rownames(meta1),
    clusterid1  = meta1$seurat_clusters
  )
  
  clusters2 <- tibble::tibble(
    cellbarcode = rownames(meta2),
    clusterid2  = meta2$seurat_clusters
  )
  
  # no renaming needed now
  clusters_merged <- dplyr::full_join(clusters1, clusters2, by = "cellbarcode")
  
  indices1 <- split(seq_len(nrow(clusters_merged)), clusters_merged$clusterid1)
  indices2 <- split(seq_len(nrow(clusters_merged)), clusters_merged$clusterid2)
  
  jaccard_list <- list()
  for (i in names(indices1)) {
    for (j in names(indices2)) {
      set1 <- indices1[[i]]
      set2 <- indices2[[j]]
      similarity <- jaccard_similarity(set1, set2)
      jaccard_list[[length(jaccard_list) + 1]] <- list(
        name1 = i,
        name2 = j,
        jaccard_similarity = similarity
      )
    }
  }
  
  jaccard_df <- do.call(rbind, lapply(jaccard_list, as.data.frame))
  jaccard_df <- type.convert(jaccard_df, as.is = TRUE)
  
  # IMPORTANT: make axes discrete (so cluster IDs show as tick labels)
  jaccard_df$name1 <- factor(jaccard_df$name1,
                             levels = sort(unique(jaccard_df$name1)))
  jaccard_df$name2 <- factor(jaccard_df$name2,
                             levels = sort(unique(jaccard_df$name2)))
  
  p <- ggplot2::ggplot(
    data = jaccard_df,
    ggplot2::aes(x = name1, y = name2, fill = jaccard_similarity)
  ) +
    ggplot2::geom_tile(color = "black") +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "firebrick",
      breaks = seq(0, 1, 0.2)
    ) +
    ggplot2::geom_point(
      data = subset(jaccard_df, jaccard_similarity >= threshold),
      ggplot2::aes(x = name1, y = name2),
      shape = 8,       # star-like marker
      size  = 1.5,     # adjust to taste
      color = "black"
    ) +
    ggplot2::scale_x_discrete(name = name1) +
    ggplot2::scale_y_discrete(name = name2) +
    ggplot2::labs(
      fill  = "Jaccard Similarity"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      axis.title.x = ggplot2::element_text(size = 14),
      axis.title.y = ggplot2::element_text(size = 14)
    )
  
  p
}

#' Compute max Jaccard similarity per reference cluster across multiple methods
#'
#' For a given reference clustering (in \code{meta1}) and a list of alternative
#' clusterings (\code{meta_list}), this function computes, for each reference
#' cluster and each method, the maximum Jaccard similarity with any cluster
#' from that method. The result is returned as a long-format data frame/tibble.
#'
#' @param meta1 A data frame (e.g. \code{seurat_obj@meta.data}) containing
#'   metadata for the reference clustering. Row names must be cell barcodes
#'   and there must be a \code{seurat_clusters} column.
#' @param meta_list A list of data frames, each analogous to \code{meta1} but
#'   corresponding to a different clustering method or parameter setting.
#'   Each element must have row names as cell barcodes and a
#'   \code{seurat_clusters} column.
#' @param method_names A character vector giving the name of each clustering
#'   method, with length equal to \code{length(meta_list)}. These names are
#'   used to label the output in the \code{method} column.
#'
#' @return A tibble (or data frame) with one row per reference cluster and
#'   method, containing:
#'   \itemize{
#'     \item \code{clusterid}: reference cluster ID (from \code{meta1}).
#'     \item \code{max_jaccard}: maximum Jaccard similarity between that
#'       reference cluster and any cluster in the corresponding method.
#'     \item \code{method}: method name, from \code{method_names}.
#'   }
#'
#' @details
#' Cell barcodes are matched between \code{meta1} and each element of
#' \code{meta_list} via a full join on row names treated as \code{cellbarcode}.
#' Jaccard similarity between clusters is computed using
#' \code{\link{jaccard_similarity}} on the sets of row indices belonging to
#' each cluster pair.
#'
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr full_join group_by summarise rename bind_rows
#'
#' @examples
#' \dontrun{
#' res <- MakeMaxJaccardInterMethodDataframe_Multi(
#'   meta1        = ref_seurat@meta.data,
#'   meta_list    = list(methodA@meta.data, methodB@meta.data),
#'   method_names = c("MethodA", "MethodB")
#' )
#' head(res)
#' }
MakeMaxJaccardInterMethodDataframe_Multi <- function(meta1, meta_list, method_names) {
  if (length(meta_list) != length(method_names)) {
    stop("meta_list and method_names must have the same length.")
  }
  
  # Precompute cluster assignments for the reference
  clusters1 <- tibble::tibble(
    cellbarcode = rownames(meta1),
    clusterid1  = meta1$seurat_clusters
  )
  
  # Helper: compute max Jaccard for one meta2
  compute_max_jaccard_one <- function(meta2) {
    clusters2 <- tibble::tibble(
      cellbarcode = rownames(meta2),
      clusterid2  = meta2$seurat_clusters
    )
    
    clusters_merged <- dplyr::full_join(clusters1, clusters2, by = "cellbarcode")
    
    indices1 <- split(seq_len(nrow(clusters_merged)), clusters_merged$clusterid1)
    indices2 <- split(seq_len(nrow(clusters_merged)), clusters_merged$clusterid2)
    
    jaccard_list <- list()
    for (i in names(indices1)) {
      for (j in names(indices2)) {
        set1 <- indices1[[i]]
        set2 <- indices2[[j]]
        similarity <- jaccard_similarity(set1, set2)
        jaccard_list[[length(jaccard_list) + 1]] <- list(
          name1 = i,
          name2 = j,
          jaccard_similarity = similarity
        )
      }
    }
    
    jaccard_df <- do.call(rbind, lapply(jaccard_list, as.data.frame))
    jaccard_df <- type.convert(jaccard_df, as.is = TRUE)
    
    jaccard_tibble <- tibble::as_tibble(jaccard_df)
    
    # One row per reference cluster with its max Jaccard to this meta2
    result <- jaccard_tibble %>%
      dplyr::group_by(name1) %>%
      dplyr::summarise(
        max_jaccard = max(jaccard_similarity),
        .groups = "drop"
      ) %>%
      dplyr::rename(clusterid = name1)
    
    result
  }
  
  # Loop over meta_list, compute, and add method column
  res_list <- vector("list", length(meta_list))
  for (k in seq_along(meta_list)) {
    res_k <- compute_max_jaccard_one(meta_list[[k]])
    res_k$method <- method_names[k]
    res_list[[k]] <- res_k
  }
  
  # Bind all methods together
  dplyr::bind_rows(res_list)
}

#' @title
#' Multi-panel boxplot of cluster stability across multiple datasets
#'
#' @description
#' Creates a faceted boxplot of cluster-wise stability (e.g. max Jaccard values)
#' for an arbitrary number of datasets, with one panel per dataset.
#'
#' @param ... One or more data frames containing per-cluster stability
#'   summaries. Each must contain at least columns \code{clusterid} and
#'   \code{max_jaccard}.
#' @param panel_labels Optional character vector giving the labels to use for
#'   the facets corresponding to the input data frames. If \code{NULL}
#'   (default), labels are set to \code{"Dataset 1"}, \code{"Dataset 2"}, etc.
#'   The length of \code{panel_labels} must match the number of data frames
#'   supplied in \code{...}.
#' @param nrow Integer; number of rows in the facet layout (default: \code{3}).
#'
#' @return A \code{ggplot2} object showing boxplots of \code{max_jaccard}
#'   by \code{clusterid}, faceted by dataset/panel. A horizontal dashed red line
#'   is drawn at \code{y = 0.75} as a visual stability threshold.
#'
#' @details
#' All input data frames are combined into a single long-format data frame
#' with an additional column \code{panel} indicating the dataset of origin.
#' Cluster IDs are treated as factors, and facets are ordered according to
#' \code{panel_labels}. This function is intended for visual comparison of
#' cluster stability distributions across multiple datasets or methods.
#'
#' @importFrom dplyr bind_rows mutate
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_hline facet_wrap labs
#'   theme_bw theme element_rect element_text
#'
#' @examples
#' \dontrun{
#' p <- MultiPanelStabilityByClusterPlot(
#'   stab_ds1, stab_ds2, stab_ds3,
#'   panel_labels = c("Dataset A", "Dataset B", "Dataset C"),
#'   nrow = 1
#' )
#' print(p)
#' }
MultiPanelStabilityByClusterPlot <- function(...,
                                             panel_labels = NULL,
                                             nrow = 3,title="Stability by cluster") {
  # Collect input data frames into a list
  df_list <- list(...)
  n_panels <- length(df_list)
  if (n_panels == 0) {
    stop("At least one data frame must be supplied.")
  }
  
  # Default labels if none supplied
  if (is.null(panel_labels)) {
    panel_labels <- paste("Dataset", seq_len(n_panels))
  }
  if (length(panel_labels) != n_panels) {
    stop("panel_labels must be NULL or a character vector with length equal to the number of data frames supplied.")
  }
  
  # Add an identifier column to each data frame and combine
  df_combined <- bind_rows(
    lapply(seq_along(df_list), function(i) {
      df <- df_list[[i]]
      df$panel <- panel_labels[i]
      df
    })
  )
  
  # Make sure clusterid is treated as discrete within each panel
  df_combined <- df_combined %>%
    mutate(
      clusterid = as.factor(clusterid),
      panel = factor(panel, levels = panel_labels)  # control facet order
    )
  
  # Build the plot
  p <- ggplot(df_combined, aes(x = clusterid, y = max_jaccard)) +
    geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white") +
    geom_hline(yintercept = 0.75, color = "red", linetype = "dashed") +
    facet_wrap(~ panel, scales = "free_x", nrow = nrow) +
    labs(
      x = "Cluster ID",
      y = expression(paste("Cluster ", stability["max jaccard"])),
      title = ""
    ) +
    theme(plot.title = element_text(hjust = 0)) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "gray", color = NA),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1)
    )
  
  return(p)
}



#' @title
#' Compare inter-method vs intra-method cluster stability (downsampling-based)
#'
#' @description
#' Computes the number of "stable" reference clusters between two clustering
#' solutions (inter-method stability) and compares this to downsampling-based
#' intra-method stability for each method, visualized as a violin/boxplot
#' with a horizontal reference line.
#'
#' @param meta1 A data frame (e.g. \code{seurat_obj@meta.data}) containing
#'   metadata for the first clustering. Row names must be cell barcodes and
#'   there must be a \code{seurat_clusters} column.
#' @param meta2 A data frame containing metadata for the second clustering.
#'   Row names must be cell barcodes and there must be a
#'   \code{seurat_clusters} column.
#' @param downsamples1 A data frame summarizing downsampling-based stability for
#'   the first method. Must contain at least the columns:
#'   \code{downsample_number}, \code{clusterid}, and \code{max_jaccard}.
#' @param downsamples2 A data frame summarizing downsampling-based stability for
#'   the second method, with the same structure as \code{downsamples1}.
#' @param threshold Numeric value in \code{[0, 1]} used as a Jaccard similarity
#'   cutoff to define a "stable" cluster both for inter-method comparison and
#'   within-method downsampling summaries.
#' @param label1 Character string; label to use on the x-axis for the first
#'   method in the final plot.
#' @param label2 Character string; label to use on the x-axis for the second
#'   method in the final plot.
#'
#' @return A \code{ggplot2} object showing, for each method, the distribution
#'   of the number of stable clusters across downsampling runs (violin + boxplot),
#'   with a dashed red line indicating the number of stable clusters between
#'   the two methods (inter-method stability).
#'
#' @details
#' Inter-method stability is computed by matching cell barcodes between
#' \code{meta1} and \code{meta2}, forming cluster-wise index sets, and
#' calculating Jaccard similarity for all pairs of clusters using
#' \code{\link{jaccard_similarity}}. The number of reference clusters (from
#' \code{meta1}) with at least one partner cluster in \code{meta2} having
#' Jaccard similarity greater than or equal to \code{threshold} is taken as
#' the inter-method stable cluster count.
#'
#' Intra-method stability for each method is derived from \code{downsamples1}
#' and \code{downsamples2} by counting, within each downsampling replicate,
#' how many clusters have \code{max_jaccard > threshold}, and then summarizing
#' these counts across downsampling runs.
#' 
#' NOTE: the Snakemake workflow that produces the downsampling files currently
#' labels the downsampling iteration in a column called \code{bootstrap_number}
#' which is why the current version of the function uses this variable when parsing
#' the files. To be changed at a later date.
#'
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr full_join group_by summarise rename bind_rows filter n_distinct
#'   left_join mutate pull
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_violin geom_hline geom_boxplot
#'   scale_x_discrete scale_y_continuous labs theme element_text
#'
#' @examples
#' \dontrun{
#' p <- MakeInterVsIntraStablePlot(
#'   meta1        = seurat_method1@meta.data,
#'   meta2        = seurat_method2@meta.data,
#'   downsamples1 = ds_res_method1,
#'   downsamples2 = ds_res_method2,
#'   threshold    = 0.5,
#'   label1       = "Method 1",
#'   label2       = "Method 2"
#' )
#' print(p)
#' }
MakeInterVsIntraStablePlot <- function(meta1, meta2,
                                       downsamples1, downsamples2,
                                       threshold, label1, label2) {
  
  clusters1 <- tibble::tibble(
    cellbarcode = rownames(meta1),
    clusterid1  = meta1$seurat_clusters
  )
  
  clusters2 <- tibble::tibble(
    cellbarcode = rownames(meta2),
    clusterid2  = meta2$seurat_clusters
  )
  
  # no renaming needed now
  clusters_merged <- dplyr::full_join(clusters1, clusters2, by = "cellbarcode")
  
  indices1 <- split(seq_len(nrow(clusters_merged)), clusters_merged$clusterid1)
  indices2 <- split(seq_len(nrow(clusters_merged)), clusters_merged$clusterid2)
  
  jaccard_list <- list()
  for (i in names(indices1)) {
    for (j in names(indices2)) {
      set1 <- indices1[[i]]
      set2 <- indices2[[j]]
      similarity <- jaccard_similarity(set1, set2)
      jaccard_list[[length(jaccard_list) + 1]] <- list(
        name1 = i,
        name2 = j,
        jaccard_similarity = similarity
      )
    }
  }
  
  jaccard_df <- do.call(rbind, lapply(jaccard_list, as.data.frame))
  jaccard_df <- type.convert(jaccard_df, as.is = TRUE)
  
  jaccard_df$name1 <- factor(jaccard_df$name1,
                             levels = sort(unique(jaccard_df$name1)))
  jaccard_df$name2 <- factor(jaccard_df$name2,
                             levels = sort(unique(jaccard_df$name2)))
  jaccard_tibble <- tibble::as_tibble(jaccard_df)
  
  stable_cluster_count <- jaccard_tibble %>%
    dplyr::filter(jaccard_similarity >= threshold) %>%
    dplyr::summarise(n = dplyr::n_distinct(name1)) %>%
    dplyr::pull(n)
  
  downsample1_summary <- downsamples1 %>%
    dplyr::group_by(bootstrap_number) %>%
    dplyr::filter(max_jaccard > threshold) %>%
    dplyr::summarise(n_stable_name1 = dplyr::n_distinct(clusterid),
                     .groups = "drop")
  
  downsample2_summary <- downsamples2 %>%
    dplyr::group_by(bootstrap_number) %>%
    dplyr::filter(max_jaccard > threshold) %>%
    dplyr::summarise(n_stable_name2 = dplyr::n_distinct(clusterid),
                     .groups = "drop")
  
  downsample_stable_merged <- downsample1_summary %>%
    dplyr::left_join(downsample2_summary, by = "bootstrap_number") %>%
    tidyr::pivot_longer(
      cols = c(n_stable_name1, n_stable_name2),
      names_to = "method",
      values_to = "n_clusters"
    ) %>%
    dplyr::mutate(method = factor(method,
                                  levels = c("n_stable_name1", "n_stable_name2")))
  
  max_y <- max(downsample_stable_merged$n_clusters, na.rm = TRUE)
  
  downsample_interstable_plot <- downsample_stable_merged %>%
    ggplot2::ggplot(ggplot2::aes(x = method, y = n_clusters)) +
    ggplot2::geom_violin(width = 0.7) +
    ggplot2::geom_hline(yintercept = stable_cluster_count,
                        color = "red", linetype = "dashed", linewidth = 1) +
    ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA) +
    ggplot2::scale_x_discrete(
      labels = c(
        n_stable_name1 = label1,
        n_stable_name2 = label2
      ),
      expand = c(0.15, 0.15)
    ) +
    ggplot2::scale_y_continuous(breaks = seq(0, max_y, by = 4)) +
    ggplot2::labs(x = "", y = "# clusters") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  return(downsample_interstable_plot)
}


#' Plot silhouette width distributions across methods
#'
#' Takes a list of metadata tables (one per method), extracts per-cell
#' silhouette widths, and produces a violin + boxplot comparison across
#' methods.
#'
#' @param meta_list A list of data frames, one per method, typically
#'   \code{seurat_obj@meta.data}. Each element must contain a
#'   \code{silhouette_width} column, and either:
#'   \itemize{
#'     \item a \code{barcode} column, or
#'     \item cell barcodes as row names (which will be converted to a
#'       \code{barcode} column).
#'   }
#' @param method_labels A character vector of method names, one per element of
#'   \code{meta_list}. Used for labeling and ordering the x-axis (and facets,
#'   if extended).
#'
#' @return A \code{ggplot2} object showing, for each method, a violin plot of
#'   the distribution of \code{silhouette_width} values with an overlaid
#'   boxplot. A horizontal red dotted line is drawn at 0.
#'
#' @details
#' For each metadata table, the function extracts \code{barcode} and
#' \code{silhouette_width}, adds a \code{method} column based on
#' \code{method_labels}, and row-binds all methods into a long-format data
#' frame. Methods are then treated as an ordered factor, and silhouette
#' distributions are visualized with \code{coord_flip()} so methods appear on
#' the y-axis.
#'
#' @importFrom purrr map2_dfr
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr select mutate
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot geom_hline
#'   coord_flip xlab ylab theme_bw
#'
#' @examples
#' \dontrun{
#' p <- SilhouetteViolinByMethodPlot(
#'   meta_list     = list(meta_method1, meta_method2, meta_method3),
#'   method_labels = c("Method 1", "Method 2", "Method 3")
#' )
#' print(p)
#' }
SilhouetteViolinByMethodPlot <- function(meta_list, method_labels) {
  if (length(meta_list) != length(method_labels)) {
    stop("meta_list and method_labels must have the same length.")
  }
  
  sil_long <- purrr::map2_dfr(
    .x = meta_list,
    .y = method_labels,
    ~ {
      meta_tbl <- .x
      method   <- .y
      
      # if barcodes are rownames, turn them into a "barcode" column
      if (!"barcode" %in% names(meta_tbl)) {
        meta_tbl <- meta_tbl %>%
          tibble::rownames_to_column(var = "barcode")
      }
      
      meta_tbl %>%
        dplyr::select(barcode, silhouette_width) %>%
        dplyr::mutate(method = method)
    }
  )
  
  sil_long <- sil_long %>%
    dplyr::mutate(method = factor(method, levels = method_labels))
  
  p <- ggplot2::ggplot(sil_long, ggplot2::aes(x = method, y = silhouette_width)) +
    ggplot2::geom_violin(trim = FALSE, fill = "grey90", color = "grey40") +
    ggplot2::geom_boxplot(width = 0.15, outlier.size = 0.3, outlier.alpha = 0.5) +
    ggplot2::geom_hline(yintercept = 0, color = "red", linetype = "dotted") +
    ggplot2::coord_flip() +
    ggplot2::xlab("Method") +
    ggplot2::ylab("Silhouette width (all barcodes)") +
    ggplot2::theme_bw()
  
  return(p)
}



#' @title
#' Plot median silhouette vs. Jaccard stability for multiple methods
#'
#' @description
#' For each clustering method, this function computes per-cluster median
#' silhouette width and per-cluster median Jaccard stability (from
#' downsampling-based summaries), then plots them as a scatter plot faceted
#' by method.
#'
#' @param meta_list A list of metadata data frames (e.g.
#'   \code{seurat_obj@meta.data}), one per method. Each element must contain
#'   columns \code{seurat_clusters} and \code{silhouette_width}.
#' @param downsamp_list A list of downsampling stability summary data frames,
#'   one per method, corresponding in order to \code{meta_list}. Each element
#'   must contain columns \code{clusterid} and \code{max_jaccard}.
#' @param method_labels A character vector of method labels, one per element
#'   of \code{meta_list} / \code{downsamp_list}. These are used for facet
#'   labels and factor levels.
#'
#' @return A \code{ggplot2} object showing, for each method, a scatter plot of
#'   per-cluster median silhouette width (x-axis) vs. per-cluster median
#'   Jaccard stability (y-axis), faceted by \code{method}.
#'
#' @details
#' For each method:
#' \itemize{
#'   \item The metadata table is grouped by \code{seurat_clusters} to compute
#'     \code{median_silhouette} (median of \code{silhouette_width} per cluster).
#'   \item The downsampling table is grouped by \code{clusterid} to compute
#'     \code{median_max_jaccard} (median of \code{max_jaccard} per cluster).
#'   \item These summaries are joined on \code{clusterid} and tagged with the
#'     corresponding \code{method} label.
#' }
#' All methods' summaries are combined into a single data frame and plotted.
#'
#' @importFrom purrr map2_dfr
#' @importFrom dplyr group_by summarise mutate inner_join
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap xlab ylab theme_bw
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' p <- ClusterStabilityVsSilhouettePlot(
#'   meta_list      = list(meta_method1, meta_method2),
#'   downsamp_list  = list(ds_method1, ds_method2),
#'   method_labels  = c("Method 1", "Method 2")
#' )
#' print(p)
#' }
ClusterStabilityVsSilhouettePlot <- function(meta_list, downsamp_list, method_labels) {
  if (length(meta_list) != length(downsamp_list) ||
      length(meta_list) != length(method_labels)) {
    stop("meta_list, downsamp_list, and method_labels must have the same length.")
  }
  
  summary_df <- purrr::map2_dfr(
    .x = seq_along(meta_list),
    .y = meta_list,
    ~ {
      i            <- .x
      meta_tbl     <- .y
      downsamp_tbl <- downsamp_list[[i]]
      method       <- method_labels[i]
      
      # column checks
      if (!all(c("seurat_clusters", "silhouette_width") %in% names(meta_tbl))) {
        stop(paste("Metadata index", i,
                   "does not have seurat_clusters and silhouette_width columns"))
      }
      if (!all(c("clusterid", "max_jaccard") %in% names(downsamp_tbl))) {
        stop(paste("Downsample index", i,
                   "does not have clusterid and max_jaccard columns"))
      }
      
      # cluster size per cluster
      size_summary <- meta_tbl %>%
        dplyr::group_by(clusterid = .data[["seurat_clusters"]]) %>%
        dplyr::summarise(
          cluster_size = dplyr::n(),
          .groups = "drop"
        ) %>%
        dplyr::mutate(clusterid = as.character(clusterid))
      
      # median silhouette per cluster
      sil_summary <- meta_tbl %>%
        dplyr::group_by(clusterid = .data[["seurat_clusters"]]) %>%
        dplyr::summarise(
          median_silhouette = median(.data[["silhouette_width"]], na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::mutate(clusterid = as.character(clusterid))
      
      # median max_jaccard per cluster
      jac_summary <- downsamp_tbl %>%
        dplyr::mutate(clusterid = as.character(clusterid)) %>%
        dplyr::group_by(clusterid) %>%
        dplyr::summarise(
          median_max_jaccard = median(max_jaccard, na.rm = TRUE),
          .groups = "drop"
        )
      
      sil_summary %>%
        dplyr::inner_join(jac_summary, by = "clusterid") %>%
        dplyr::inner_join(size_summary, by = "clusterid") %>%
        dplyr::mutate(method = method)
    }
  )
  
  summary_df <- summary_df %>%
    dplyr::mutate(
      clusterid = factor(clusterid),
      method    = factor(method, levels = method_labels)
    )
  
  ggplot(summary_df,
         aes(x = median_silhouette,
             y = median_max_jaccard,
             fill = cluster_size)) +          # map to fill, not color
    geom_point(
      shape = 21,        # filled circle with border
      color = "black",   # border color
      size  = 2.0,       # slightly larger than default
      alpha = 0.8
    ) +
    facet_wrap(~ method, scales = "fixed") +
    scale_x_continuous(limits = c(-1, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_fill_gradient(
      name = "Cluster size",
      low  = "mistyrose",   # light red
      high = "firebrick"          # strong red
    ) +
    xlab("Median silhouette width") +
    ylab(expression(paste("Cluster ", stability[Jaccard]))) +
    theme_bw()
}


### Cell type annotation

#' @title
#' Add CellTypist annotations to a Seurat object and save (with strict barcode checks)
#'
#' @param seu        A Seurat object (used to generate input for CellTypist)
#' @param labels_csv Path to CellTypist CSV file with columns:
#'                   barcode,predicted_labels,majority_voting,
#'                   prob_predicted_labels,prob_majority_voting
#' @param output_rds Path to save updated Seurat object (RDS)
#'
#' @return The updated Seurat object (invisibly)
CelltypistAnnots2SeuratMetadata <- function(seu,
                                            labels_csv,
                                            output_rds) {
  # Read CellTypist labels; 'barcode' column becomes rownames
  labels <- read.csv(
    labels_csv,
    row.names   = 1,
    check.names = FALSE
  )
  
  # Barcode sets
  seu_barcodes    <- colnames(seu)
  label_barcodes  <- rownames(labels)
  
  # Check for mismatches
  missing_in_labels <- setdiff(seu_barcodes, label_barcodes)
  extra_in_labels   <- setdiff(label_barcodes, seu_barcodes)
  
  if (length(missing_in_labels) > 0) {
    stop(
      "These Seurat barcodes are missing from labels_csv (showing up to 10): ",
      paste(head(missing_in_labels, 10), collapse = ", "),
      if (length(missing_in_labels) > 10) " ..." else ""
    )
  }
  
  if (length(extra_in_labels) > 0) {
    warning(
      "These barcodes are present in labels_csv but not in the Seurat object (showing up to 10): ",
      paste(head(extra_in_labels, 10), collapse = ", "),
      if (length(extra_in_labels) > 10) " ..." else ""
    )
  }
  
  # Align labels to Seurat barcodes
  labels_aligned <- labels[seu_barcodes, , drop = FALSE]
  
  # Add to metadata
  seu@meta.data <- cbind(
    seu@meta.data,
    labels_aligned
  )
  
  # Save updated Seurat object
  saveRDS(seu, file = output_rds)
  
  invisible(seu)
}
