args <- commandArgs(trailingOnly = TRUE)
output_path <- args[1]
marker_paths <- args[-1]

library("tidyverse")

marker_tables <- lapply(
  marker_paths,
  function(marker_path) read_csv(marker_path, show_col_types = FALSE)
)

combined_markers <- bind_rows(marker_tables)
write_csv(combined_markers, file = output_path)
