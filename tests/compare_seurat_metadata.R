args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0 || length(args) %% 2 != 0) {
  stop("Provide current/reference RDS path pairs", call. = FALSE)
}

suppressPackageStartupMessages(library("Seurat"))

numeric_tolerance <- as.numeric(Sys.getenv("SEURAT_METADATA_NUMERIC_TOLERANCE", "1e-8"))
if (is.na(numeric_tolerance) || length(numeric_tolerance) != 1 || numeric_tolerance < 0) {
  stop("SEURAT_METADATA_NUMERIC_TOLERANCE must be a non-negative number", call. = FALSE)
}

neighborhood_purity_tolerance <- as.numeric(Sys.getenv("SEURAT_METADATA_NEIGHBORHOOD_PURITY_TOLERANCE", "0.025"))
if (is.na(neighborhood_purity_tolerance) || length(neighborhood_purity_tolerance) != 1 || neighborhood_purity_tolerance < 0) {
  stop("SEURAT_METADATA_NEIGHBORHOOD_PURITY_TOLERANCE must be a non-negative number", call. = FALSE)
}

metadata_column_tolerances <- c(
  neighborhood_purity = neighborhood_purity_tolerance
)

# Every difference is collected and reported at the end, so that one run of the full test-data workflow
# (which takes hours) shows all of them, not just the first.
failures <- character(0)
failed_objects <- character(0)

fail <- function(path, message) {
  failures <<- c(failures, sprintf("%s: %s", path, message))
  failed_objects <<- union(failed_objects, path)
  invisible(FALSE)
}

compare_metadata_column <- function(current_path, column, current_values, reference_values) {
  current_na <- is.na(current_values)
  reference_na <- is.na(reference_values)
  if (!identical(current_na, reference_na)) {
    return(fail(current_path, sprintf("metadata column %s has different NA positions", column)))
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
      column_tolerance <- if (column %in% names(metadata_column_tolerances)) {
        metadata_column_tolerances[[column]]
      } else {
        numeric_tolerance
      }
      bad <- diff > column_tolerance * scale
      if (any(bad)) {
        return(fail(
          current_path,
          sprintf(
            "metadata column %s differs numerically in %d of %d cells; max abs diff %.12g (tolerance %.3g, scaled by max(|value|, 1))",
            column,
            sum(bad),
            length(bad),
            max(diff),
            column_tolerance
          )
        ))
      }
    }
    return(invisible(TRUE))
  }

  current_character <- as.character(current_values)
  reference_character <- as.character(reference_values)
  comparable <- !current_na & !reference_na
  if (!identical(current_character[comparable], reference_character[comparable])) {
    return(fail(
      current_path,
      sprintf(
        "metadata column %s differs in %d of %d cells",
        column,
        sum(current_character[comparable] != reference_character[comparable]),
        sum(comparable)
      )
    ))
  }

  invisible(TRUE)
}

compare_metadata_pair <- function(current_path, reference_path) {
  current <- readRDS(current_path)
  reference <- readRDS(reference_path)

  current_metadata <- current@meta.data
  reference_metadata <- reference@meta.data

  if (!identical(colnames(current_metadata), colnames(reference_metadata))) {
    return(fail(
      current_path,
      paste0(
        "metadata columns differ\ncurrent: ",
        paste(colnames(current_metadata), collapse = ","),
        "\nreference: ",
        paste(colnames(reference_metadata), collapse = ",")
      )
    ))
  }

  if (!identical(rownames(current_metadata), rownames(reference_metadata))) {
    return(fail(
      current_path,
      sprintf(
        "metadata barcodes differ (%d cells, reference %d; %d in common)",
        nrow(current_metadata),
        nrow(reference_metadata),
        length(intersect(rownames(current_metadata), rownames(reference_metadata)))
      )
    ))
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
  tryCatch(
    compare_metadata_pair(args[[i]], args[[i + 1]]),
    error = function(e) fail(args[[i]], paste("could not be compared:", conditionMessage(e)))
  )
}

if (length(failures) > 0) {
  message(sprintf(
    "Seurat metadata differs from the reference in %d of %d object(s):",
    length(failed_objects),
    length(args) / 2
  ))
  message(paste(failures, collapse = "\n"))
  quit(save = "no", status = 1)
}

cat(sprintf("Compared Seurat metadata for %d object(s)\n", length(args) / 2))
