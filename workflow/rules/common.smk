import os
from pathlib import Path


REQUIRED_CONFIG_KEYS = {
    "conda-channel-priority",
    "sampleTable",
    "workflow_seed",
    "emptydrop_removal_methods",
    "ambient_decon_methods",
    "doublet_removal_methods",
    "posthoc_methods",
    "min_nfeature",
    "min_ncount",
    "max_mtdna",
}
ALLOWED_EMPTYDROP_METHODS = {"tenx", "emptydrops"}
ALLOWED_DECON_METHODS = {"soupx", "cellbender_fromraw"}
ALLOWED_DOUBLET_METHODS = {"doubletfinder", "scdblfinder"}
ALLOWED_POSTHOC_METHODS = {"threshold", "mad"}


def require_non_empty_string(config_values, key, errors):
    value = config_values.get(key)
    if not isinstance(value, str) or not value.strip():
        errors.append(f"{key} must be a non-empty string")


def require_positive_int(config_values, key, errors):
    value = config_values.get(key)
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        errors.append(f"{key} must be a positive integer")


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


def validate_workflow_config(config_values):
    errors = []
    missing_keys = sorted(REQUIRED_CONFIG_KEYS - set(config_values.keys()))
    if missing_keys:
        errors.append("missing required key(s): " + ", ".join(missing_keys))

    require_non_empty_string(config_values, "sampleTable", errors)
    if config_values.get("conda-channel-priority") != "strict":
        errors.append("conda-channel-priority must be 'strict'")

    workflow_seed = config_values.get("workflow_seed")
    if isinstance(workflow_seed, bool) or not isinstance(workflow_seed, int):
        errors.append("workflow_seed must be an integer")

    emptydrop_methods = validate_method_list(config_values, "emptydrop_removal_methods", ALLOWED_EMPTYDROP_METHODS, errors)
    decon_methods = validate_method_list(config_values, "ambient_decon_methods", ALLOWED_DECON_METHODS, errors)
    doublet_methods = validate_method_list(config_values, "doublet_removal_methods", ALLOWED_DOUBLET_METHODS, errors)
    posthoc_methods = validate_method_list(config_values, "posthoc_methods", ALLOWED_POSTHOC_METHODS, errors)

    require_positive_int(config_values, "min_nfeature", errors)
    require_positive_int(config_values, "min_ncount", errors)
    max_mtdna = config_values.get("max_mtdna")
    if isinstance(max_mtdna, bool) or not isinstance(max_mtdna, (int, float)) or not 0 <= max_mtdna <= 100:
        errors.append("max_mtdna must be a number between 0 and 100")

    if "resultsDir" in config_values:
        require_non_empty_string(config_values, "resultsDir", errors)

    if errors:
        raise ValueError("Invalid workflow config: " + "; ".join(errors))

    return {
        "emptydrop_methods": emptydrop_methods,
        "decon_methods": decon_methods,
        "doublet_methods": doublet_methods,
        "posthoc_methods": posthoc_methods,
    }


REQUIRED_SAMPLE_COLUMNS = {"sampleid", "tenx_datadir"}
REQUIRED_TENX_PATHS = (
    "filtered_feature_bc_matrix",
    "raw_feature_bc_matrix",
    "raw_feature_bc_matrix.h5",
)


def resolve_tenx_datadir(path_value):
    path = Path(str(path_value)).expanduser()
    if not path.is_absolute():
        path = Path.cwd() / path
    return path


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

    if errors:
        raise ValueError(f"Invalid sample sheet {sample_table}: " + "; ".join(errors))

    validated = sampleinfo.copy()
    validated["sampleid"] = sample_ids.astype(str)
    validated["tenx_datadir"] = resolved_data_dirs
    return validated


def sample_tenx_dir(wildcards):
    return sampleinfo.loc[sampleinfo["sampleid"] == wildcards.sample, "tenx_datadir"].values[0]


def tenx_filtered_matrix_input(wildcards):
    return os.path.join(sample_tenx_dir(wildcards), "filtered_feature_bc_matrix")


def tenx_raw_matrix_input(wildcards):
    return os.path.join(sample_tenx_dir(wildcards), "raw_feature_bc_matrix")


def tenx_raw_h5_input(wildcards):
    return os.path.join(sample_tenx_dir(wildcards), "raw_feature_bc_matrix.h5")


def marker_chunk_inputs(wildcards):
    manifest_dir = checkpoints.marker_manifest.get(prefix=wildcards.prefix).output.manifest
    cluster_ids = glob_wildcards(str(Path(manifest_dir) / "{cluster}.txt")).cluster
    return expand(
        f"{RESULTS_DIR}/{{prefix}}_markergenes_cluster{{cluster}}.csv",
        prefix=wildcards.prefix,
        cluster=cluster_ids,
    )
