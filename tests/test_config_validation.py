import copy
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml


CONFIG_FILE = Path("config/config.yaml")
TEST_SAMPLE_SHEET = Path("testdata/samplesheet_test.tsv")
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
ALLOWED_WORKFLOW_MODES = {"preprocess", "preprocess_and_downsample", "downsample_only"}


def load_default_config(repo_root):
    with (repo_root / CONFIG_FILE).open() as handle:
        return yaml.safe_load(handle)


def write_config(path, config):
    with path.open("w") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)
    return path


def results_dir_path(tmp_path):
    return tmp_path / "results"


def run_snakemake(repo_root, config_file, results_dir, *extra_args):
    snakemake = shutil.which("snakemake")
    assert snakemake is not None, "snakemake is not available on PATH"

    cmd = [
        snakemake,
        *extra_args,
        "--profile",
        "none",
        "--workflow-profile",
        "none",
        "--snakefile",
        "workflow/Snakefile",
        "--configfile",
        str(config_file),
    ]
    if results_dir is not None:
        cmd.extend(["--config", f"resultsDir={results_dir}"])
    return subprocess.run(
        cmd,
        cwd=repo_root,
        text=True,
        capture_output=True,
        check=False,
        timeout=120,
    )


def combined_output(result):
    return result.stdout + result.stderr


def baseline_test_config(repo_root):
    config = copy.deepcopy(load_default_config(repo_root))
    config["sampleTable"] = TEST_SAMPLE_SHEET.as_posix()
    return config


def assert_no_duplicate_values(config, key):
    values = config[key]
    assert len(values) == len(set(values)), f"{key} contains duplicates"


def test_default_config_has_required_keys_and_valid_values():
    repo_root = Path(__file__).resolve().parents[1]
    config = load_default_config(repo_root)

    assert REQUIRED_CONFIG_KEYS <= set(config)
    assert config["conda-channel-priority"] == "strict"
    assert isinstance(config["sampleTable"], str) and config["sampleTable"]
    assert config.get("workflow_mode") in ALLOWED_WORKFLOW_MODES
    assert isinstance(config["workflow_seed"], int) and not isinstance(config["workflow_seed"], bool)
    assert isinstance(config["nDownsampleReplicates"], int) and config["nDownsampleReplicates"] > 0
    assert isinstance(config["downsampleRate"], (int, float)) and 0 < config["downsampleRate"] <= 1
    assert isinstance(config["downsampleSeuratObjectDir"], str) and config["downsampleSeuratObjectDir"]
    assert isinstance(config["downsampleResultsDir"], str) and config["downsampleResultsDir"]

    for key, allowed_values in [
        ("emptydrop_removal_methods", ALLOWED_EMPTYDROP_METHODS),
        ("ambient_decon_methods", ALLOWED_DECON_METHODS),
        ("doublet_removal_methods", ALLOWED_DOUBLET_METHODS),
        ("posthoc_methods", ALLOWED_POSTHOC_METHODS),
    ]:
        assert isinstance(config[key], list) and config[key]
        assert set(config[key]) <= allowed_values
        assert_no_duplicate_values(config, key)

    assert isinstance(config["min_nfeature"], int) and config["min_nfeature"] > 0
    assert isinstance(config["min_ncount"], int) and config["min_ncount"] > 0
    assert isinstance(config["max_mtdna"], (int, float)) and 0 <= config["max_mtdna"] <= 100

    assert config.get("preflight_mode", "skip") in {"off", "warn", "skip", "error"}
    if "preflight_min_cells" in config:
        assert isinstance(config["preflight_min_cells"], int) and config["preflight_min_cells"] > 0


@pytest.mark.parametrize(
    "mutate, expected_message, override_results_dir",
    [
        (lambda cfg: cfg.pop("posthoc_methods"), "missing required key(s): posthoc_methods", True),
        (lambda cfg: cfg.update({"posthoc_methods": ["threshold", "bad_method"]}), "posthoc_methods contains invalid value(s): bad_method", True),
        (lambda cfg: cfg.update({"doublet_removal_methods": ["scdblfinder", "scdblfinder"]}), "doublet_removal_methods contains duplicate value(s): scdblfinder", True),
        (lambda cfg: cfg.update({"ambient_decon_methods": []}), "ambient_decon_methods must be a non-empty list", True),
        (lambda cfg: cfg.update({"workflow_seed": "12345"}), "workflow_seed must be an integer", True),
        (lambda cfg: cfg.update({"min_nfeature": 0}), "min_nfeature must be a positive integer", True),
        (lambda cfg: cfg.update({"max_mtdna": 101}), "max_mtdna must be a number between 0 and 100", True),
        (lambda cfg: cfg.update({"resultsDir": ""}), "resultsDir must be a non-empty string", False),
        (lambda cfg: cfg.update({"workflow_mode": "bad_mode"}), "workflow_mode must be one of", True),
        (lambda cfg: cfg.update({"workflow_mode": "downsample_only", "downsampleRate": 1.5}), "downsampleRate must be > 0 and <= 1", True),
        (lambda cfg: cfg.update({"preflight_min_cells": 0}), "preflight_min_cells must be a positive integer", True),
        (lambda cfg: cfg.update({"preflight_mode": "bogus"}), "preflight_mode must be one of", True),
        (lambda cfg: cfg.update({"cellbender_learning_rate": 0}), "cellbender_learning_rate must be a positive number", True),
        (lambda cfg: cfg.update({"cellbender_learning_rate": "fast"}), "cellbender_learning_rate must be a positive number", True),
        (lambda cfg: cfg.update({"cellbender_adaptive_rerun": "yes"}), "cellbender_adaptive_rerun must be a boolean", True),
    ],
)
def test_invalid_config_fails_early_with_clear_message(tmp_path, mutate, expected_message, override_results_dir):
    repo_root = Path(__file__).resolve().parents[1]
    config = baseline_test_config(repo_root)
    mutate(config)
    config_file = write_config(tmp_path / "invalid_config.yaml", config)

    results_dir = results_dir_path(tmp_path) if override_results_dir else None
    result = run_snakemake(repo_root, config_file, results_dir, "-np")
    output = combined_output(result)

    assert result.returncode != 0
    assert "Invalid workflow config" in output
    assert expected_message in output


def run_variant_dry_run(tmp_path, **overrides):
    repo_root = Path(__file__).resolve().parents[1]
    config = baseline_test_config(repo_root)
    config.update(overrides)
    config_file = write_config(tmp_path / "variant_config.yaml", config)
    results_dir = results_dir_path(tmp_path)
    result = run_snakemake(repo_root, config_file, results_dir, "-np")
    output = combined_output(result)
    assert result.returncode == 0, output
    return output, results_dir


def test_posthoc_method_variant_changes_rule_all_targets(tmp_path):
    output, results_dir = run_variant_dry_run(tmp_path, posthoc_methods=["threshold"])

    assert str(results_dir / "posthocfilter" / "seurat_posthocfilt_threshold_scdblfinder_cellbender_fromraw_test_markergenes.csv") in output
    assert str(results_dir / "posthocfilter" / "seurat_posthocfilt_threshold_doubletfinder_soupx_tenx_test_markergenes.csv") in output
    assert "seurat_posthocfilt_mad" not in output


def test_doublet_method_variant_changes_rule_all_targets(tmp_path):
    output, results_dir = run_variant_dry_run(tmp_path, doublet_removal_methods=["scdblfinder"])

    assert str(results_dir / "scdblfinder" / "seurat_scdblfinder_soupx_tenx_test_markergenes.csv") in output
    assert str(results_dir / "scdblfinder" / "seurat_scdblfinder_cellbender_fromraw_test_markergenes.csv") in output
    assert "seurat_doubletfinder" not in output
    assert "doubletfinder_installed.txt" not in output


def test_decon_method_variant_can_exclude_cellbender_targets(tmp_path):
    output, results_dir = run_variant_dry_run(tmp_path, ambient_decon_methods=["soupx"])

    assert str(results_dir / "soupx" / "seurat_soupx_tenx_test_markergenes.csv") in output
    assert str(results_dir / "soupx" / "seurat_soupx_emptydrops_test_markergenes.csv") in output
    assert "cellbender" not in output


def test_emptydrop_method_variant_excludes_emptydrops_specific_targets(tmp_path):
    output, results_dir = run_variant_dry_run(tmp_path, emptydrop_removal_methods=["tenx"])

    assert str(results_dir / "seurat_filtered" / "filtered_seurat_tenx_test_markergenes.csv") in output
    assert str(results_dir / "soupx" / "seurat_soupx_tenx_test_markergenes.csv") in output
    assert "emptydrops" not in output


def test_downsample_only_mode_builds_dag_from_external_seurat_objects(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    config = {
        "conda-channel-priority": "strict",
        "workflow_mode": "downsample_only",
        "downsampleSeuratObjectDir": "testdata/downsampling/seurat_objects",
        "downsampleResultsDir": (tmp_path / "downsampling").as_posix(),
        "nDownsampleReplicates": 2,
        "downsampleRate": 0.5,
        "workflowSeed": 12345,
    }
    config_file = write_config(tmp_path / "downsample_only.yaml", config)

    result = run_snakemake(repo_root, config_file, None, "-np")
    output = combined_output(result)

    assert result.returncode == 0, output
    assert "downsample_clusters" in output
    assert str(tmp_path / "downsampling" / "filtered_seurat_tenx_test_clusterdownsampling.tsv") in output
    assert "tenx2seuratrds" not in output


def test_preprocess_and_downsample_mode_adds_downsample_targets(tmp_path):
    output, results_dir = run_variant_dry_run(
        tmp_path,
        workflow_mode="preprocess_and_downsample",
        downsampleTargets=["filtered_seurat_tenx_test"],
        downsampleResultsDir=(tmp_path / "downsampling").as_posix(),
    )

    assert "tenx2seuratrds" in output
    assert "downsample_clusters" in output
    assert str(tmp_path / "downsampling" / "filtered_seurat_tenx_test_clusterdownsampling.tsv") in output
