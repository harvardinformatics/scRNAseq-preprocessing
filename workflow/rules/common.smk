import gzip
import os
import re
from pathlib import Path


BASE_REQUIRED_CONFIG_KEYS = {
    "conda-channel-priority",
}
PREPROCESS_REQUIRED_CONFIG_KEYS = {
    "sampleTable",
    "workflow_seed",
    "emptydrop_removal_methods",
    "ambient_decon_methods",
    "doublet_removal_methods",
    "posthoc_methods",
    "min_nfeature",
    "min_ncount",
}
ALLOWED_WORKFLOW_MODES = {"preprocess", "preprocess_and_downsample", "downsample_only"}
ALLOWED_EMPTYDROP_METHODS = {"tenx", "emptydrops"}
ALLOWED_DECON_METHODS = {"soupx", "cellbender_fromraw", "decontx_raw"}
ALLOWED_DOUBLET_METHODS = {"doubletfinder", "scdblfinder"}
ALLOWED_POSTHOC_METHODS = {"threshold", "mad"}
ALLOWED_PREFLIGHT_MODES = {"off", "warn", "skip", "error"}
DEFAULT_PREFLIGHT_MIN_CELLS = 100
DEFAULT_PREFLIGHT_MODE = "skip"
DEFAULT_MIN_RAW_TO_CELL_RATIO = 2
# CellBender's own default learning rate; the adaptive re-run uses exactly half of this.
DEFAULT_CELLBENDER_LEARNING_RATE = 0.0001
DEFAULT_CELLBENDER_ADAPTIVE_RERUN = True
# The cellbendercells object (CellBender's cell calls on raw counts) is opt-in.
DEFAULT_CELLBENDER_CELL_CALLS = False
# Cell metadata columns every SCTransform call may regress out (its vars.to.regress); none by default.
ALLOWED_VARS_TO_REGRESS = {"percent.mt", "nCount_RNA", "nFeature_RNA"}
DEFAULT_VARS_TO_REGRESS = []
# Without a sample-sheet mt_genes list, percent.mt counts the features whose Seurat row names match the R
# pattern "(?i)^mt-" (workflow/scripts/regression_utils.R); this is the same test in Python.
MT_NAME_PATTERN = re.compile(r"mt-", re.IGNORECASE)


def require_non_empty_string(config_values, key, errors):
    value = config_values.get(key)
    if not isinstance(value, str) or not value.strip():
        errors.append(f"{key} must be a non-empty string")


def require_positive_int(config_values, key, errors):
    value = config_values.get(key)
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        errors.append(f"{key} must be a positive integer")


def require_optional_positive_int(config_values, key, errors):
    if key in config_values:
        require_positive_int(config_values, key, errors)


def require_optional_positive_number(config_values, key, errors):
    if key not in config_values:
        return
    value = config_values.get(key)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or value <= 0:
        errors.append(f"{key} must be a positive number")


def require_optional_bool(config_values, key, errors):
    if key in config_values and not isinstance(config_values.get(key), bool):
        errors.append(f"{key} must be a boolean (true or false)")


def validate_method_list(config_values, key, allowed_values, errors):
    value = config_values.get(key)
    if not isinstance(value, list) or not value:
        errors.append(f"{key} must be a non-empty list")
        return []

    normalized = []
    for item in value:
        if not isinstance(item, str) or not item.strip():
            errors.append(f"{key} contains a non-string or empty value")
            continue
        normalized.append(item)

    duplicates = sorted({item for item in normalized if normalized.count(item) > 1})
    if duplicates:
        errors.append(f"{key} contains duplicate value(s): " + ", ".join(duplicates))

    invalid = sorted(set(normalized) - allowed_values)
    if invalid:
        errors.append(
            f"{key} contains invalid value(s): "
            + ", ".join(invalid)
            + "; allowed values are: "
            + ", ".join(sorted(allowed_values))
        )

    return normalized


def validate_vars_to_regress(config_values, errors):
    """The optional vars_to_regress list (default []): metadata columns SCTransform regresses out."""
    if "vars_to_regress" not in config_values:
        return list(DEFAULT_VARS_TO_REGRESS)
    allowed = ", ".join(sorted(ALLOWED_VARS_TO_REGRESS))
    value = config_values["vars_to_regress"]
    if not isinstance(value, list):
        errors.append(f"vars_to_regress must be a list (default []) of any of: {allowed}")
        return []
    normalized = []
    for item in value:
        if not isinstance(item, str) or not item.strip():
            errors.append("vars_to_regress contains a non-string or empty value")
            continue
        normalized.append(item)
    duplicates = sorted({item for item in normalized if normalized.count(item) > 1})
    if duplicates:
        errors.append("vars_to_regress contains duplicate value(s): " + ", ".join(duplicates))
    invalid = sorted(set(normalized) - ALLOWED_VARS_TO_REGRESS)
    if invalid:
        errors.append(f"vars_to_regress contains invalid value(s): {', '.join(invalid)}; allowed values are: {allowed}")
    return normalized


def validated_vars_to_regress(config_values):
    """validate_vars_to_regress() for entry points that do not run validate_workflow_config()."""
    errors = []
    value = validate_vars_to_regress(config_values, errors)
    if errors:
        raise ValueError("Invalid workflow config: " + "; ".join(errors))
    return value


def validate_downsample_targets(config_values, errors):
    if "downsampleTargets" not in config_values:
        return

    value = config_values["downsampleTargets"]
    if value == "all":
        return
    if not isinstance(value, list) or not value:
        errors.append("downsampleTargets must be 'all' or a non-empty list")
        return
    for item in value:
        if not isinstance(item, str) or not item.strip():
            errors.append("downsampleTargets contains a non-string or empty value")
            return

    duplicates = sorted({item for item in value if value.count(item) > 1})
    if duplicates:
        errors.append("downsampleTargets contains duplicate value(s): " + ", ".join(duplicates))


def validate_excluded_samples(config_values, errors):
    """Normalize the optional user-controlled excluded_samples list (default empty).

    Sample IDs listed here are dropped from the DAG. Exclusion is explicit and visible so
    nothing is skipped automatically (e.g. a re-sequenced library reusing a sample ID is
    processed normally unless the user deliberately lists that ID here)."""
    value = config_values.get("excluded_samples", [])
    if value is None:
        return []
    if not isinstance(value, list):
        errors.append("excluded_samples must be a list of sample id strings")
        return []
    normalized = []
    for item in value:
        if not isinstance(item, str) or not item.strip():
            errors.append("excluded_samples contains a non-string or empty value")
            continue
        normalized.append(item.strip())
    return sorted(set(normalized))


def validate_workflow_config(config_values):
    errors = []
    workflow_mode = config_values.get("workflow_mode", "preprocess")
    if workflow_mode not in ALLOWED_WORKFLOW_MODES:
        errors.append(
            "workflow_mode must be one of: " + ", ".join(sorted(ALLOWED_WORKFLOW_MODES))
        )
        workflow_mode = "preprocess"

    required_keys = set(BASE_REQUIRED_CONFIG_KEYS)
    if workflow_mode in {"preprocess", "preprocess_and_downsample"}:
        required_keys.update(PREPROCESS_REQUIRED_CONFIG_KEYS)

    missing_keys = sorted(required_keys - set(config_values.keys()))
    if missing_keys:
        errors.append("missing required key(s): " + ", ".join(missing_keys))

    if config_values.get("conda-channel-priority") != "strict":
        errors.append("conda-channel-priority must be 'strict'")

    emptydrop_methods = []
    decon_methods = []
    doublet_methods = []
    posthoc_methods = []
    if workflow_mode in {"preprocess", "preprocess_and_downsample"}:
        require_non_empty_string(config_values, "sampleTable", errors)

        workflow_seed = config_values.get("workflow_seed")
        if isinstance(workflow_seed, bool) or not isinstance(workflow_seed, int):
            errors.append("workflow_seed must be an integer")

        emptydrop_methods = validate_method_list(config_values, "emptydrop_removal_methods", ALLOWED_EMPTYDROP_METHODS, errors)
        decon_methods = validate_method_list(config_values, "ambient_decon_methods", ALLOWED_DECON_METHODS, errors)
        doublet_methods = validate_method_list(config_values, "doublet_removal_methods", ALLOWED_DOUBLET_METHODS, errors)
        posthoc_methods = validate_method_list(config_values, "posthoc_methods", ALLOWED_POSTHOC_METHODS, errors)

        require_positive_int(config_values, "min_nfeature", errors)
        require_positive_int(config_values, "min_ncount", errors)
        if "max_mtdna" in config_values:
            errors.append("max_mtdna is no longer supported: the pipeline does not filter cells on % mitochondrial "
                          "reads (mitochondrial genes are recognisable by name in few species). percent.mt is "
                          "kept as cell metadata; to regress it out in SCTransform, add percent.mt to "
                          "vars_to_regress")

        require_optional_positive_int(config_values, "preflight_min_cells", errors)
        if config_values.get("preflight_mode", DEFAULT_PREFLIGHT_MODE) not in ALLOWED_PREFLIGHT_MODES:
            errors.append(
                "preflight_mode must be one of: " + ", ".join(sorted(ALLOWED_PREFLIGHT_MODES))
            )

        require_optional_positive_number(config_values, "cellbender_learning_rate", errors)
        require_optional_bool(config_values, "cellbender_adaptive_rerun", errors)
        require_optional_bool(config_values, "cellbender_cell_calls", errors)
        if config_values.get("cellbender_cell_calls") is True and "cellbender_fromraw" not in decon_methods:
            errors.append("cellbender_cell_calls needs cellbender_fromraw in ambient_decon_methods: "
                          "its CellBender run supplies the cell calls")

    if workflow_mode in {"preprocess_and_downsample", "downsample_only"}:
        if "downsampleSeuratObjectDir" in config_values:
            require_non_empty_string(config_values, "downsampleSeuratObjectDir", errors)
        if "downsampleResultsDir" in config_values:
            require_non_empty_string(config_values, "downsampleResultsDir", errors)
        require_optional_positive_int(config_values, "nDownsampleReplicates", errors)
        if "workflowSeed" in config_values:
            seed = config_values["workflowSeed"]
            if isinstance(seed, bool) or not isinstance(seed, int):
                errors.append("workflowSeed must be an integer")
        if "downsampleRate" in config_values:
            rate = config_values["downsampleRate"]
            if isinstance(rate, bool) or not isinstance(rate, (int, float)) or not 0 < rate <= 1:
                errors.append("downsampleRate must be > 0 and <= 1")
        validate_downsample_targets(config_values, errors)

    if "resultsDir" in config_values:
        require_non_empty_string(config_values, "resultsDir", errors)

    excluded_samples = validate_excluded_samples(config_values, errors)
    # SCTransform runs in preprocessing and in every downsampling replicate, so this applies in every mode.
    vars_to_regress = validate_vars_to_regress(config_values, errors)

    if errors:
        raise ValueError("Invalid workflow config: " + "; ".join(errors))

    return {
        "workflow_mode": workflow_mode,
        "emptydrop_methods": emptydrop_methods,
        "decon_methods": decon_methods,
        "doublet_methods": doublet_methods,
        "posthoc_methods": posthoc_methods,
        "excluded_samples": excluded_samples,
        "vars_to_regress": vars_to_regress,
        "preflight_min_cells": config_values.get("preflight_min_cells", DEFAULT_PREFLIGHT_MIN_CELLS),
        "preflight_mode": config_values.get("preflight_mode", DEFAULT_PREFLIGHT_MODE),
        "cellbender_learning_rate": config_values.get(
            "cellbender_learning_rate", DEFAULT_CELLBENDER_LEARNING_RATE
        ),
        "cellbender_adaptive_rerun": config_values.get(
            "cellbender_adaptive_rerun", DEFAULT_CELLBENDER_ADAPTIVE_RERUN
        ),
        "cellbender_cell_calls": config_values.get(
            "cellbender_cell_calls", DEFAULT_CELLBENDER_CELL_CALLS
        ),
    }


REQUIRED_SAMPLE_COLUMNS = {"sampleid", "tenx_datadir"}
# Optional: a file of the sample's mitochondrial gene IDs; blank (or no column) means the "(?i)^mt-" pattern.
MT_GENES_COLUMN = "mt_genes"
REQUIRED_TENX_PATHS = (
    "filtered_feature_bc_matrix",
    "raw_feature_bc_matrix",
    "raw_feature_bc_matrix.h5",
)


def resolve_sample_sheet_path(path_value):
    """A sample-sheet path: absolute, or relative to the working directory (~ expanded)."""
    path = Path(str(path_value)).expanduser()
    if not path.is_absolute():
        path = Path.cwd() / path
    return path


def resolve_tenx_datadir(path_value):
    return resolve_sample_sheet_path(path_value)


def read_mt_gene_ids(path):
    """The gene IDs in an mt_genes file, in file order without duplicates: one per line, surrounding spaces,
    tabs and carriage returns trimmed, blank lines and lines starting with # skipped (as regression_utils.R)."""
    ids = {}
    with open(path) as handle:
        for line in handle:
            line = line.strip(" \t\r\n")
            if line and not line.startswith("#"):
                ids.setdefault(line, None)
    return list(ids)


def first_few(items, n=5):
    shown = ", ".join(items[:n])
    return f"{shown}, ... ({len(items) - n} more)" if len(items) > n else shown


def tenx_features_file(matrix_dir):
    """The features file Seurat's Read10X reads: genes.tsv (Cell Ranger < 3) if present, else features.tsv.gz."""
    genes = Path(matrix_dir) / "genes.tsv"
    return genes if genes.exists() else Path(matrix_dir) / "features.tsv.gz"


def read_tenx_features(matrix_dir):
    """(path, ids, names): columns 1 and 2 of a 10x matrix directory's features file, one per feature."""
    path = tenx_features_file(matrix_dir)
    opener = gzip.open if path.name.endswith(".gz") else open
    ids = []
    names = []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t")
            ids.append(fields[0].strip(" \t\r\n"))
            names.append(fields[1] if len(fields) > 1 else "")
    return path, ids, names


def seurat_row_name_start(feature_id, name):
    """How the Seurat row name of a feature starts: Read10X takes column 2 (column 1 where it is empty or NA)
    and make.unique only appends suffixes; CreateSeuratObject turns "_" and "|" into "-"."""
    if name in ("", "NA"):
        name = feature_id
    return name.replace("_", "-").replace("|", "-")


def mt_gene_checks(sampleinfo, vars_to_regress):
    """Parse-time mitochondrial-gene checks, reading only each sample's filtered features file. Returns
    (errors, warnings):
      errors   -- an mt_genes list none of whose IDs is in column 1 of the features; or, with percent.mt in
                  vars_to_regress, a sample without a list whose feature names have none matching "(?i)^mt-".
      warnings -- listed IDs that are not features (they are skipped).
    The R scripts check the same at run time, on the objects they build."""
    errors = []
    warnings = []
    regress_mt = "percent.mt" in vars_to_regress
    for sample_id, data_dir, mt_genes in zip(
        sampleinfo["sampleid"], sampleinfo["tenx_datadir"], sampleinfo[MT_GENES_COLUMN]
    ):
        if not mt_genes and not regress_mt:
            continue
        matrix_dir = Path(data_dir) / "filtered_feature_bc_matrix"
        try:
            features_path, feature_ids, names = read_tenx_features(matrix_dir)
        except (OSError, UnicodeDecodeError) as exc:
            errors.append(f"{sample_id}: cannot read the features of {matrix_dir}: {exc}")
            continue
        if mt_genes:
            ids = read_mt_gene_ids(mt_genes)
            known = set(feature_ids)
            missing = [gene_id for gene_id in ids if gene_id not in known]
            if len(missing) == len(ids):
                errors.append(
                    f"{sample_id}: none of the {len(ids)} gene IDs in its mt_genes file {mt_genes} is in "
                    f"column 1 of {features_path} (e.g. {first_few(missing)})"
                )
            elif missing:
                warnings.append(
                    f"{sample_id}: {len(missing)} of {len(ids)} gene IDs in its mt_genes file {mt_genes} are "
                    f"not in column 1 of {features_path} and are skipped: {first_few(missing)}"
                )
        elif not any(
            MT_NAME_PATTERN.match(seurat_row_name_start(feature_id, name))
            for feature_id, name in zip(feature_ids, names)
        ):
            errors.append(
                f"{sample_id}: percent.mt is in vars_to_regress, but no feature name in {features_path} "
                "matches \"(?i)^mt-\", so percent.mt would be 0 in every cell. List the sample's mitochondrial "
                "gene IDs in a file and give its path in the sample sheet's mt_genes column, or remove "
                "percent.mt from vars_to_regress"
            )
    return errors, warnings


def validate_sample_sheet(sampleinfo, sample_table):
    errors = []
    missing_columns = sorted(REQUIRED_SAMPLE_COLUMNS - set(sampleinfo.columns))
    if missing_columns:
        errors.append("missing required column(s): " + ", ".join(missing_columns))
        raise ValueError(f"Invalid sample sheet {sample_table}: " + "; ".join(errors))

    if sampleinfo.empty:
        errors.append("sample sheet has no rows")

    sample_ids = sampleinfo["sampleid"].astype("string")
    missing_sample_ids = sampleinfo.index[sample_ids.isna() | (sample_ids.str.strip() == "")].tolist()
    if missing_sample_ids:
        errors.append("sampleid is missing or empty on row(s): " + ", ".join(str(i + 2) for i in missing_sample_ids))

    duplicate_sample_ids = sorted(sample_ids[sample_ids.duplicated(keep=False)].dropna().unique())
    if duplicate_sample_ids:
        errors.append("duplicate sampleid value(s): " + ", ".join(duplicate_sample_ids))

    data_dirs = sampleinfo["tenx_datadir"].astype("string")
    missing_data_dirs = sampleinfo.index[data_dirs.isna() | (data_dirs.str.strip() == "")].tolist()
    if missing_data_dirs:
        errors.append("tenx_datadir is missing or empty on row(s): " + ", ".join(str(i + 2) for i in missing_data_dirs))

    if errors:
        raise ValueError(f"Invalid sample sheet {sample_table}: " + "; ".join(errors))

    resolved_data_dirs = []
    for row_number, (sample_id, data_dir_value) in enumerate(zip(sample_ids, data_dirs), start=2):
        data_dir = resolve_tenx_datadir(data_dir_value)
        resolved_data_dirs.append(str(data_dir))
        if not data_dir.exists():
            errors.append(f"sample {sample_id} row {row_number}: tenx_datadir does not exist: {data_dir}")
            continue
        for required_path in REQUIRED_TENX_PATHS:
            candidate = data_dir / required_path
            if not candidate.exists():
                errors.append(f"sample {sample_id} row {row_number}: missing {required_path}: {candidate}")

    # Optional mt_genes column: blank means the "(?i)^mt-" pattern; a path (resolved like tenx_datadir) must
    # be a file listing at least one gene ID.
    if MT_GENES_COLUMN in sampleinfo.columns:
        mt_values = sampleinfo[MT_GENES_COLUMN].astype("string")
    else:
        mt_values = pd.Series([pd.NA] * len(sampleinfo), dtype="string")
    resolved_mt_genes = []
    for row_number, (sample_id, mt_value) in enumerate(zip(sample_ids, mt_values), start=2):
        if pd.isna(mt_value) or not mt_value.strip():
            resolved_mt_genes.append("")
            continue
        mt_path = resolve_sample_sheet_path(mt_value.strip())
        resolved_mt_genes.append(str(mt_path))
        if not mt_path.is_file():
            errors.append(f"sample {sample_id} row {row_number}: mt_genes file does not exist: {mt_path}")
            continue
        try:
            mt_ids = read_mt_gene_ids(mt_path)
        except (OSError, UnicodeDecodeError) as exc:
            errors.append(f"sample {sample_id} row {row_number}: mt_genes file is not readable text: {mt_path} ({exc})")
            continue
        if not mt_ids:
            errors.append(
                f"sample {sample_id} row {row_number}: mt_genes file lists no gene IDs "
                f"(blank lines and lines starting with # are ignored): {mt_path}"
            )

    if errors:
        raise ValueError(f"Invalid sample sheet {sample_table}: " + "; ".join(errors))

    validated = sampleinfo.copy()
    validated["sampleid"] = sample_ids.astype(str)
    validated["tenx_datadir"] = resolved_data_dirs
    validated[MT_GENES_COLUMN] = resolved_mt_genes
    return validated


def count_filtered_cells(tenx_datadir):
    """Number of called cells in a sample's CellRanger filtered matrix.

    Reads only the barcode list (barcodes.tsv[.gz]) line count -- no matrix is loaded -- so it
    is cheap enough to run for every sample at parse time (the preflight check below)."""
    matrix_dir = Path(tenx_datadir) / "filtered_feature_bc_matrix"
    for name in ("barcodes.tsv.gz", "barcodes.tsv"):
        barcodes = matrix_dir / name
        if barcodes.exists():
            opener = gzip.open if name.endswith(".gz") else open
            with opener(barcodes, "rt") as handle:
                return sum(1 for line in handle if line.strip())
    raise FileNotFoundError(f"preflight cell count: no barcodes.tsv[.gz] under {matrix_dir}")


def preflight_min_cells_check(sampleinfo, min_cells):
    """Per-sample CellRanger-filtered cell counts vs the preflight minimum.

    Returns (counts, below): counts maps each sampleid to its number of called cells in
    filtered_feature_bc_matrix; below is the sorted list of sampleids with fewer than
    min_cells. Only the CellRanger filtered count is knowable before the run (emptydrops and
    cellbender_fromraw call cells at runtime), so this gates on that count as the whole-sample
    viability signal; require_min_cells_for_pca remains the runtime backstop for the arms."""
    counts = {
        str(sample_id): count_filtered_cells(data_dir)
        for sample_id, data_dir in zip(sampleinfo["sampleid"], sampleinfo["tenx_datadir"])
    }
    below = sorted(sample_id for sample_id, n in counts.items() if n < min_cells)
    return counts, below


def read_mtx_dims(matrix_dir):
    """(n_features, n_barcodes, nnz) from a 10x MatrixMarket matrix.mtx[.gz] header.

    Reads only the header (comment lines plus the single dims line), so it never loads the
    matrix -- cheap enough for every sample at parse time."""
    for name in ("matrix.mtx.gz", "matrix.mtx"):
        mtx = Path(matrix_dir) / name
        if mtx.exists():
            opener = gzip.open if name.endswith(".gz") else open
            with opener(mtx, "rt") as handle:
                for line in handle:
                    if line.startswith("%"):
                        continue
                    parts = line.split()
                    if len(parts) < 3:
                        raise ValueError(f"malformed MatrixMarket dims line: {line.strip()!r}")
                    return int(parts[0]), int(parts[1]), int(parts[2])
            raise ValueError(f"no dimension line in {mtx}")
    raise FileNotFoundError(f"no matrix.mtx[.gz] under {matrix_dir}")


def cellbender_preflight_checks(sampleinfo, min_raw_to_cell_ratio=DEFAULT_MIN_RAW_TO_CELL_RATIO):
    """Raw-matrix sanity for the cellbender_fromraw arm, at parse time (reads only MTX headers).

    Returns (errors, warnings). These are input-integrity checks, deliberately NOT governed by
    preflight_mode:
      errors   -- a broken or wrong raw input to fix: the raw matrix.mtx is unreadable,
                  empty/degenerate (zero features/barcodes/nonzeros), or has fewer barcodes than
                  the sample's filtered matrix (raw must be a superset of filtered). Hard-stops.
      warnings -- the raw matrix has fewer than min_raw_to_cell_ratio x the called cells in
                  droplets: too few empty droplets for CellBender's ambient estimate, often a
                  sign the raw path actually points at a filtered matrix. Non-blocking.
    Only meaningful when cellbender_fromraw is a configured decon method; the caller gates on that."""
    errors = []
    warnings = []
    for sample_id, data_dir in zip(sampleinfo["sampleid"], sampleinfo["tenx_datadir"]):
        raw_dir = Path(data_dir) / "raw_feature_bc_matrix"
        try:
            n_features, n_barcodes, nnz = read_mtx_dims(raw_dir)
        except Exception as exc:
            errors.append(f"{sample_id}: raw matrix unreadable ({raw_dir}): {exc}")
            continue
        if n_features <= 0 or n_barcodes <= 0 or nnz <= 0:
            errors.append(
                f"{sample_id}: raw matrix is empty/degenerate "
                f"(features={n_features}, barcodes={n_barcodes}, nonzeros={nnz})"
            )
            continue
        filtered_cells = count_filtered_cells(data_dir)
        if n_barcodes < filtered_cells:
            errors.append(
                f"{sample_id}: raw matrix has fewer barcodes ({n_barcodes}) than the filtered "
                f"matrix ({filtered_cells}); the raw and filtered inputs look mismatched"
            )
        elif n_barcodes < min_raw_to_cell_ratio * filtered_cells:
            warnings.append(
                f"{sample_id}: {n_barcodes} raw droplets vs {filtered_cells} called cells "
                f"(< {min_raw_to_cell_ratio}x); few empty droplets for CellBender's ambient "
                "estimate - check the raw path is a true unfiltered matrix"
            )
    return errors, warnings


def sample_tenx_dir(wildcards):
    return sampleinfo.loc[sampleinfo["sampleid"] == wildcards.sample, "tenx_datadir"].values[0]


def tenx_filtered_matrix_input(wildcards):
    return os.path.join(sample_tenx_dir(wildcards), "filtered_feature_bc_matrix")


def tenx_raw_matrix_input(wildcards):
    return os.path.join(sample_tenx_dir(wildcards), "raw_feature_bc_matrix")


def tenx_raw_h5_input(wildcards):
    return os.path.join(sample_tenx_dir(wildcards), "raw_feature_bc_matrix.h5")


def sample_mt_genes(wildcards):
    """The sample's mt_genes file (absolute path), or "" for the "(?i)^mt-" pattern: the scripts' --mt-genes."""
    return str(sampleinfo.loc[sampleinfo["sampleid"] == wildcards.sample, MT_GENES_COLUMN].values[0])


def sample_mt_genes_input(wildcards):
    """The sample's mt_genes file as a rule input, so that editing the list reruns the steps that read it."""
    mt_genes = sample_mt_genes(wildcards)
    return [mt_genes] if mt_genes else []


def decontx_cells_matrix_input(wildcards):
    """Called-cell (filtered) matrix for the decontx_raw arm, selected by empty_method.

    The decontx_raw arm decontaminates the same called-cell sets as soupx:
    the CellRanger filtered matrix (tenx) or the emptydrops filtered matrix. Mirrors the filtered
    inputs of soupx.smk (tenx) and soupx_emptydrops.smk (emptydrops)."""
    if wildcards.empty_method == "emptydrops":
        return f"{RESULTS_DIR}/emptydrops/{wildcards.sample}_emptydrops_filtered_matrix"
    return tenx_filtered_matrix_input(wildcards)


def decontx_seurat_base_input(wildcards):
    """Base Seurat object whose clusters decontX reuses as its cell populations (z), by empty_method."""
    if wildcards.empty_method == "emptydrops":
        return f"{RESULTS_DIR}/emptydrops/filtered_seurat_emptydrops_{wildcards.sample}.rds"
    return f"{RESULTS_DIR}/seurat_filtered/filtered_seurat_tenx_{wildcards.sample}.rds"


def marker_chunk_inputs(wildcards):
    manifest_dir = checkpoints.marker_manifest.get(prefix=wildcards.prefix).output.manifest
    cluster_ids = glob_wildcards(str(Path(manifest_dir) / "{cluster}.txt")).cluster
    return expand(
        f"{RESULTS_DIR}/{{prefix}}_markergenes_cluster{{cluster}}.csv",
        prefix=wildcards.prefix,
        cluster=cluster_ids,
    )


def marker_targets(prefixes):
    return [f"{RESULTS_DIR}/{prefix}_markergenes.csv" for prefix in prefixes]


def rds_targets(prefixes):
    return [f"{RESULTS_DIR}/{prefix}.rds" for prefix in prefixes]


def select_downsample_inputs(inputs_by_target):
    requested = config.get("downsampleTargets", "all")
    if requested == "all" or requested == ["all"]:
        return dict(inputs_by_target)

    missing = sorted(set(requested) - set(inputs_by_target))
    if missing:
        raise ValueError(
            "downsampleTargets contains unknown target(s): "
            + ", ".join(missing)
            + "; available targets are: "
            + ", ".join(sorted(inputs_by_target))
        )
    return {target: inputs_by_target[target] for target in requested}


def downsample_inputs_from_preprocess_outputs():
    inputs = {}
    for path in PREPROCESS_SEURAT_TARGETS:
        target = Path(path).stem
        if target in inputs and inputs[target] != path:
            raise ValueError(f"duplicate downsample target name: {target}")
        inputs[target] = path
    return inputs


def downsample_inputs_from_external_dir():
    targets = sorted(
        glob_wildcards(f"{DOWNSAMPLE_SEURAT_OBJECT_DIR}/{{downsample_target}}.rds").downsample_target
    )
    return {
        target: f"{DOWNSAMPLE_SEURAT_OBJECT_DIR}/{target}.rds"
        for target in targets
    }


def finalize_run(results_dir, all_targets, samplesheet, enabled):
    """Post-run quarantine + completeness verification that drives the process exit code.

    Runs from the Snakefile's onsuccess/onerror handlers. It (1) quarantines low-quality
    samples and (2) verifies that every required target exists, failing the run only when a
    required output is missing for a sample that was NOT flagged low quality. sys.exit() in a
    handler deterministically sets Snakemake's exit code either way, so the runner batch job's
    COMPLETED/FAILED state reflects the true outcome even when Snakemake would exit 0 on a
    terminal failure.

    `enabled` should be True only for non-local (SLURM) runs: the verification exists to
    correct that executor's unreliable exit code, whereas local runs (tests, ad-hoc builds)
    have reliable exit codes and may intentionally build only a subset of targets.
    """
    if not enabled:
        return

    import subprocess
    import sys as _sys
    from pathlib import Path as _Path

    _Path(results_dir).mkdir(parents=True, exist_ok=True)
    targets_file = _Path(results_dir) / ".run_targets.txt"
    targets_file.write_text("\n".join(all_targets) + "\n")

    subprocess.run([
        _sys.executable, "workflow/scripts/quarantine_low_quality_samples.py",
        "--results-dir", results_dir, "--samplesheet", samplesheet,
    ])
    result = subprocess.run([
        _sys.executable, "workflow/scripts/verify_run_complete.py",
        "--results-dir", results_dir, "--samplesheet", samplesheet,
        "--targets-file", str(targets_file),
    ])
    _sys.exit(result.returncode)
