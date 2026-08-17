args <- commandArgs(trailingOnly = TRUE)
output_path <- args[1]
marker_paths <- args[-1]

library("tidyverse")

marker_tables <- lapply(
  marker_paths,
  function(marker_path) read_csv(marker_path, show_col_types = FALSE)
)

# Clusters with < 3 cells produce header-only (0-row) files. read_csv infers their
# empty columns as logical, which would clash with the character/double types from
# populated files under bind_rows. Drop empties before binding; if every cluster was
# skipped, fall back to a single header-only table so the output still has the schema.
non_empty <- marker_tables[vapply(marker_tables, nrow, integer(1)) > 0]
combined_markers <- if (length(non_empty) > 0) {
  bind_rows(non_empty)
} else {
  marker_tables[[1]]
}
write_csv(combined_markers, file = output_path)
