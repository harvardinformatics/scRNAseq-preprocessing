import csv
import math
import os
import shutil
import subprocess
from pathlib import Path


TEST_SEURAT_DIR = Path("testdata/downsampling/seurat_objects")
TEST_RESULTS_DIR = Path("testdata/downsampling/results")
REFERENCE_ROOT = Path("tests/downsampling/reference_outputs")
REFERENCE_RESULTS_DIR = Path("testdata/results")
TEST_SAMPLES = [
    "filtered_seurat_emptydrops_test",
    "filtered_seurat_tenx_test",
]
TEST_CONFIG = [
    "workflow_mode=downsample_only",
    f"downsampleSeuratObjectDir={TEST_SEURAT_DIR.as_posix()}",
    f"downsampleResultsDir={TEST_RESULTS_DIR.as_posix()}",
    "nDownsampleReplicates=2",
    "downsampleRate=0.5",
    "workflowSeed=12345",
]
EXPECTED_COLUMNS = ["clusterid", "max_jaccard", "bootstrap_number"]
EXPECTED_BOOTSTRAPS = {1, 2}


def repo_root():
    return Path(__file__).resolve().parents[2]


def snakemake_executable():
    snakemake = shutil.which("snakemake")
    assert snakemake is not None, "snakemake is not available on PATH"
    return snakemake


def run_command(cmd, cwd, timeout=300):
    return subprocess.run(
        cmd,
        cwd=cwd,
        text=True,
        capture_output=True,
        check=False,
        timeout=timeout,
    )


def combined_output(result):
    return result.stdout + result.stderr


def base_snakemake_cmd():
    return [
        snakemake_executable(),
        "--snakefile",
        "workflow/Snakefile",
        "--configfile",
        "config/config.yaml",
        "--config",
        *TEST_CONFIG,
        "--profile",
        "none",
        "--workflow-profile",
        "none",
        "--executor",
        "local",
    ]


def workflow_run_cmd():
    return base_snakemake_cmd() + [
        "--cores",
        "1",
        "--jobs",
        "1",
        "--latency-wait",
        "30",
        "--rerun-incomplete",
        "--use-conda",
    ]


def expected_output_for_sample(sample, results_dir=TEST_RESULTS_DIR):
    return Path(results_dir) / f"{sample}_clusterdownsampling.tsv"


def reference_output_for_sample(sample):
    return REFERENCE_ROOT / REFERENCE_RESULTS_DIR / f"{sample}_clusterdownsampling.tsv"


def read_tsv(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        return reader.fieldnames or [], rows


def assert_stability_table(path):
    columns, rows = read_tsv(path)
    assert rows, f"{path}: downsample output is empty"
    assert columns == EXPECTED_COLUMNS

    bootstraps = {int(row["bootstrap_number"]) for row in rows}
    assert bootstraps == EXPECTED_BOOTSTRAPS
    for row in rows:
        assert row["clusterid"] != ""
        max_jaccard = float(row["max_jaccard"])
        assert 0 <= max_jaccard <= 1


def rows_by_bootstrap_cluster(path, rows):
    indexed = {}
    for row in rows:
        key = (int(row["bootstrap_number"]), row["clusterid"])
        assert key not in indexed, f"{path}: duplicate bootstrap/cluster row: {key}"
        indexed[key] = row
    return indexed


def assert_stability_table_matches_reference(current_path, reference_path):
    current_columns, current_rows = read_tsv(current_path)
    reference_columns, reference_rows = read_tsv(reference_path)
    assert current_columns == reference_columns == EXPECTED_COLUMNS

    current_index = rows_by_bootstrap_cluster(current_path, current_rows)
    reference_index = rows_by_bootstrap_cluster(reference_path, reference_rows)
    assert set(current_index) == set(reference_index), (
        f"{current_path}: bootstrap/cluster rows differ from reference"
    )

    rel_tol = float(os.getenv("DOWNSAMPLE_REFERENCE_REL_TOLERANCE", "1e-6"))
    abs_tol = float(os.getenv("DOWNSAMPLE_REFERENCE_ABS_TOLERANCE", "1e-8"))
    for key, reference_row in reference_index.items():
        current_row = current_index[key]
        assert current_row["bootstrap_number"] == reference_row["bootstrap_number"]
        assert current_row["clusterid"] == reference_row["clusterid"]
        current_jaccard = float(current_row["max_jaccard"])
        reference_jaccard = float(reference_row["max_jaccard"])
        assert math.isclose(
            current_jaccard,
            reference_jaccard,
            rel_tol=rel_tol,
            abs_tol=abs_tol,
        ), (
            f"{current_path}: max_jaccard differs for bootstrap/cluster {key}: "
            f"{current_jaccard} != {reference_jaccard}"
        )


def assert_workflow_outputs_match_references(root):
    for sample in TEST_SAMPLES:
        current = root / expected_output_for_sample(sample)
        reference = root / reference_output_for_sample(sample)
        assert reference.exists(), f"missing reference output: {reference}"
        assert_stability_table_matches_reference(current, reference)
