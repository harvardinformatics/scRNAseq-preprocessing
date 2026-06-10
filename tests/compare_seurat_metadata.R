args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || length(args) %% 2 != 0) {
  stop("Provide current/reference RDS path pairs", call. = FALSE)
}

suppressPackageStartupMessages(library("Seurat"))

numeric_tolerance <- as.numeric(Sys.getenv("SEURAT_METADATA_NUMERIC_TOLERANCE", "1e-8"))
if (is.na(numeric_tolerance) || length(numeric_tolerance) != 1) {
  stop("SEURAT_METADATA_NUMERIC_TOLERANCE must be numeric", call. = FALSE)
}

fail <- function(path, message) {
  stop(sprintf("%s: %s", path, message), call. = FALSE)
}

compare_metadata_column <- function(current_path, column, current_values, reference_values) {
  current_na <- is.na(current_values)
  reference_na <- is.na(reference_values)
  if (!identical(current_na, reference_na)) {
    fail(current_path, sprintf("metadata column %s has different NA positions", column))
  }

  both_numeric <- is.numeric(current_values) && is.numeric(reference_values)
  both_integer <- is.integer(current_values) && is.integer(reference_values)
  if (both_numeric || both_integer) {
    current_numeric <- as.numeric(current_values)
    reference_numeric <- as.numeric(reference_values)
    comparable <- !current_na & !reference_na
    if (any(comparable)) {
      diff <- abs(current_numeric[comparable] - reference_numeric[comparable])
      scale <- pmax(abs(current_numeric[comparable]), abs(reference_numeric[comparable]), 1)
      bad <- diff > numeric_tolerance * scale
      if (any(bad)) {
        fail(
          current_path,
          sprintf(
            "metadata column %s differs numerically; max abs diff %.12g",
            column,
            max(diff)
          )
        )
      }
    }
    return(invisible(TRUE))
  }

  current_character <- as.character(current_values)
  reference_character <- as.character(reference_values)
  comparable <- !current_na & !reference_na
  if (!identical(current_character[comparable], reference_character[comparable])) {
    fail(current_path, sprintf("metadata column %s differs", column))
  }

  invisible(TRUE)
}

compare_metadata_pair <- function(current_path, reference_path) {
  current <- readRDS(current_path)
  reference <- readRDS(reference_path)

  current_metadata <- current@meta.data
  reference_metadata <- reference@meta.data

  if (!identical(colnames(current_metadata), colnames(reference_metadata))) {
    fail(
      current_path,
      paste0(
        "metadata columns differ\ncurrent: ",
        paste(colnames(current_metadata), collapse = ","),
        "\nreference: ",
        paste(colnames(reference_metadata), collapse = ",")
      )
    )
  }

  if (!identical(rownames(current_metadata), rownames(reference_metadata))) {
    fail(current_path, "metadata barcodes differ")
  }

  for (column in colnames(current_metadata)) {
    compare_metadata_column(
      current_path,
      column,
      current_metadata[[column]],
      reference_metadata[[column]]
    )
  }

  invisible(TRUE)
}

for (i in seq(1, length(args), by = 2)) {
  compare_metadata_pair(args[[i]], args[[i + 1]])
}

cat(sprintf("Compared Seurat metadata for %d object(s)\n", length(args) / 2))
