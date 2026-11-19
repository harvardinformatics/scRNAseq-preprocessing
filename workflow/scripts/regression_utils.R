# percent.mt and SCTransform's regression, shared by every script that computes percent.mt or calls
# SCTransform (config vars_to_regress, sample-sheet column mt_genes; see README, "Mitochondrial genes and
# regression"). The options follow a script's positional arguments as --name=value:
#   --vars-to-regress=A,B  cell metadata columns SCTransform regresses out (its vars.to.regress): any of
#                          percent.mt, nCount_RNA, nFeature_RNA. Empty or absent: none (NULL), the default.
#   --mt-genes=FILE        the sample's mitochondrial gene IDs, one per line (blank lines and lines starting
#                          with # are ignored), matched against column 1 of its features. Empty or absent:
#                          percent.mt counts the features whose names match MT_PATTERN.
#   --sample=NAME          the sample, named in messages.
# With every option empty or absent, a script computes percent.mt with MT_PATTERN and calls SCTransform with
# vars.to.regress = NULL (its default), as before these options existed.
#
# percent.mt is, in every saved object, the share of each cell's RNA counts (nCount_RNA) in its mitochondrial
# features. Every script that saves an object computes it with add_percent_mt(), on the RNA assay whatever the
# object's default assay, so it is never computed on an SCT assay (the doublet-removal scripts' inputs default
# to theirs) nor carried over from an upstream object that computed it there.
#
# A list is mapped to row names once, by the script that builds the object from the sample's features: a
# listed ID selects the rows whose column-1 ID it is, and the row names are taken from the object itself, so
# they are whatever Seurat named those rows (Read10X's make.unique of column 2, the "_" and "|" that
# CreateSeuratObject turns into "-", emptyDrops' make.unique of the Symbol, CellBender's h5 names). The
# names are stored in the object (misc slot MT_FEATURES_MISC) and every later step that starts from that
# object reads them from there, because a saved object no longer has its feature IDs.

MT_PATTERN <- "(?i)^mt-"
ALLOWED_VARS_TO_REGRESS <- c("percent.mt", "nCount_RNA", "nFeature_RNA")
MT_FEATURES_MISC <- "percent_mt_features"
SCRIPT_OPTIONS <- c("vars-to-regress", "mt-genes", "sample")

# Split the command line into positional arguments and the options above.
parse_script_args <- function(args) {
  is_option <- startsWith(args, "--")
  options <- args[is_option]
  bad <- options[!grepl("=", options, fixed = TRUE)]
  if (length(bad) > 0) {
    stop("options take the form --name=value: ", paste(bad, collapse = " "), call. = FALSE)
  }
  keys <- sub("=.*$", "", substring(options, 3))
  values <- trimws(sub("^[^=]*=", "", options))
  unknown <- setdiff(keys, SCRIPT_OPTIONS)
  if (length(unknown) > 0) {
    stop("unknown option(s): ", paste0("--", unknown, collapse = ", "), "; known: ",
         paste0("--", SCRIPT_OPTIONS, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(keys)) {
    stop("option(s) given more than once: ", paste0("--", unique(keys[duplicated(keys)]), collapse = ", "),
         call. = FALSE)
  }
  value_of <- function(key) {
    value <- values[keys == key]
    if (length(value) == 0 || !nzchar(value)) NULL else value
  }
  list(
    positional = args[!is_option],
    vars_to_regress = parse_vars_to_regress(value_of("vars-to-regress")),
    mt_genes = value_of("mt-genes"),
    sample = value_of("sample")
  )
}

# NULL for no regression, else the validated variable names.
parse_vars_to_regress <- function(value) {
  if (is.null(value)) {
    return(NULL)
  }
  vars <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  if (length(vars) == 0 || !all(nzchar(vars))) {
    stop("--vars-to-regress has an empty name: ", value, call. = FALSE)
  }
  unknown <- setdiff(vars, ALLOWED_VARS_TO_REGRESS)
  if (length(unknown) > 0) {
    stop("--vars-to-regress: unknown variable(s) ", paste(unknown, collapse = ", "), "; allowed: ",
         paste(ALLOWED_VARS_TO_REGRESS, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(vars)) {
    stop("--vars-to-regress names a variable more than once: ", value, call. = FALSE)
  }
  vars
}

message_prefix <- function(opts, context) {
  if (is.null(opts$sample)) context else sprintf("%s (sample %s)", context, opts$sample)
}

first_few <- function(x, n = 5) {
  shown <- paste(utils::head(x, n), collapse = ", ")
  if (length(x) > n) sprintf("%s, ... (%d more)", shown, length(x) - n) else shown
}

# The IDs in an mt_genes file: one per line, surrounding spaces, tabs and carriage returns trimmed, blank
# lines and lines starting with # skipped, duplicates dropped.
read_mt_gene_ids <- function(path) {
  if (!file.exists(path) || dir.exists(path)) {
    stop("mt_genes file does not exist: ", path, call. = FALSE)
  }
  lines <- trimws(readLines(path, warn = FALSE))
  ids <- unique(lines[nzchar(lines) & !startsWith(lines, "#")])
  if (length(ids) == 0) {
    stop("mt_genes file lists no gene IDs (blank lines and lines starting with # are ignored): ", path,
         call. = FALSE)
  }
  ids
}

# Column 1 (the feature IDs) of the features file Seurat::Read10X reads from a 10x matrix directory -
# genes.tsv (Cell Ranger < 3) if present, else features.tsv.gz - in row order.
read_feature_ids <- function(matrix_dir) {
  path <- file.path(matrix_dir, "genes.tsv")
  if (!file.exists(path)) {
    path <- file.path(matrix_dir, "features.tsv.gz")
  }
  if (!file.exists(path)) {
    stop("no genes.tsv or features.tsv.gz in ", matrix_dir, call. = FALSE)
  }
  con <- gzfile(path, "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  lines <- lines[nzchar(lines)]
  structure(trimws(sub("\t.*$", "", lines)), source = path)
}

# The feature IDs of a CellBender (10x v3 layout) h5, in the row order of scCustomize::Read_CellBender_h5_Mat().
read_h5_feature_ids <- function(h5) {
  ids <- tryCatch(
    as.character(rhdf5::h5read(h5, "matrix/features/id")),
    error = function(e) {
      stop("cannot read the feature IDs (matrix/features/id) of ", h5, ": ", conditionMessage(e), call. = FALSE)
    }
  )
  structure(trimws(ids), source = paste0(h5, " (matrix/features/id)"))
}

# With an mt_genes list, the files listing the same features in the same order is what lets one file's IDs
# name another matrix's rows; stop if they do not.
require_same_features <- function(row_names, other_row_names, opts, context, what) {
  if (!identical(row_names, other_row_names)) {
    stop(sprintf("%s: cannot map the mt_genes IDs: %s", message_prefix(opts, context), what), call. = FALSE)
  }
  invisible(TRUE)
}

# The row names of `object` whose feature ID (feature_ids, one per row in row order) is in the sample's
# mt_genes list. Reports the listed IDs that are not features and stops if none is.
map_mt_genes <- function(opts, feature_ids, object, context) {
  prefix <- message_prefix(opts, context)
  ids <- read_mt_gene_ids(opts$mt_genes)
  source <- attr(feature_ids, "source")
  if (length(feature_ids) != nrow(object)) {
    stop(sprintf("%s: cannot map the mt_genes IDs: %s lists %d features, the object has %d rows",
                 prefix, source, length(feature_ids), nrow(object)), call. = FALSE)
  }
  missing <- setdiff(ids, feature_ids)
  if (length(missing) == length(ids)) {
    stop(sprintf("%s: none of the %d gene IDs in the mt_genes file %s is in column 1 of %s (e.g. %s)",
                 prefix, length(ids), opts$mt_genes, source, first_few(missing)), call. = FALSE)
  }
  if (length(missing) > 0) {
    message(sprintf("%s: %d of %d mt_genes IDs are not in column 1 of %s and are skipped: %s",
                    prefix, length(missing), length(ids), source, first_few(missing)))
  }
  features <- rownames(object)[feature_ids %in% ids]
  message(sprintf("%s: %d mt_genes ID(s) from %s name %d feature(s): %s",
                  prefix, length(ids) - length(missing), opts$mt_genes, length(features), first_few(features)))
  features
}

# The mitochondrial row names an upstream step mapped from the sample's mt_genes list, or NULL when it used
# MT_PATTERN.
stored_mt_features <- function(object) {
  object@misc[[MT_FEATURES_MISC]]
}

# The mitochondrial features of the RNA assay: the listed ones it has, or those matching MT_PATTERN.
rna_mt_features <- function(object, mt_features) {
  rows <- rownames(object[["RNA"]])
  if (is.null(mt_features)) grep(pattern = MT_PATTERN, x = rows, value = TRUE) else intersect(mt_features, rows)
}

# Regressing percent.mt needs at least one mitochondrial feature; otherwise percent.mt is 0 in every cell
# and the regression silently does nothing.
require_mt_features <- function(object, mt_features, opts, context) {
  if (length(rna_mt_features(object, mt_features)) > 0) {
    return(invisible(TRUE))
  }
  prefix <- message_prefix(opts, context)
  if (is.null(mt_features)) {
    stop(sprintf(paste0(
      "%s: percent.mt is in vars_to_regress, but no feature of the RNA assay matches \"%s\", so percent.mt ",
      "is 0 in every cell. List the sample's mitochondrial gene IDs in a file and give its path in the ",
      "sample sheet's mt_genes column, or remove percent.mt from vars_to_regress."),
      prefix, MT_PATTERN), call. = FALSE)
  }
  stop(sprintf(paste0(
    "%s: percent.mt is in vars_to_regress, but none of the %d features mapped from the sample's mt_genes ",
    "list is in the RNA assay (e.g. %s), so percent.mt is 0 in every cell. Check the mt_genes list, or ",
    "remove percent.mt from vars_to_regress."),
    prefix, length(mt_features), first_few(mt_features)), call. = FALSE)
}

# percent.mt with PercentageFeatureSet() on the RNA assay, whatever the object's default assay: the share of
# each cell's RNA counts (nCount_RNA) in the listed features (mt_features, row names) or, when there is no list
# (NULL), in the features whose names match MT_PATTERN. Regressing percent.mt out needs a mitochondrial feature.
# A list is stored in the object for the steps that start from it. These are the workflow's only calls of
# PercentageFeatureSet(), each with assay = "RNA" (tests/test_vars_to_regress.py checks both).
add_percent_mt <- function(object, mt_features, opts, context) {
  regressed <- "percent.mt" %in% opts$vars_to_regress
  if (regressed) {
    require_mt_features(object, mt_features, opts, context)
  }
  if (is.null(mt_features)) {
    if (regressed) {
      message(sprintf("%s: percent.mt from the %d feature(s) matching \"%s\" (RNA assay)", message_prefix(opts, context),
                      length(rna_mt_features(object, NULL)), MT_PATTERN))
    }
    object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = MT_PATTERN, assay = "RNA")
  } else {
    used <- rna_mt_features(object, mt_features)
    message(sprintf("%s: percent.mt from %d of the sample's %d mt_genes feature(s) (RNA assay)",
                    message_prefix(opts, context), length(used), length(mt_features)))
    object[["percent.mt"]] <- PercentageFeatureSet(object, features = used, assay = "RNA")
    object@misc[[MT_FEATURES_MISC]] <- mt_features
  }
  object
}

# Before SCTransform(vars.to.regress = opts$vars_to_regress): every variable must be a cell metadata column.
check_regression_vars <- function(object, opts, context) {
  vars <- opts$vars_to_regress
  if (is.null(vars)) {
    return(invisible(NULL))
  }
  missing <- setdiff(vars, colnames(object[[]]))
  if (length(missing) > 0) {
    stop(sprintf("%s: cannot regress %s out: not in the cell metadata", message_prefix(opts, context),
                 paste(missing, collapse = ", ")), call. = FALSE)
  }
  message(sprintf("%s: SCTransform regresses out %s", message_prefix(opts, context), paste(vars, collapse = ", ")))
  invisible(vars)
}
