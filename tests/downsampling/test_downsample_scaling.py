import time

import pytest

from utils import (
    EXPECTED_COLUMNS,
    combined_output,
    read_tsv,
    repo_root,
    run_command,
    snakemake_executable,
)


SMALL_SAMPLE = "filtered_seurat_tenx_test"
SMALL_NCELLS = 300
MEDIUM_SAMPLE = "filtered_seurat_tenx_test_medium"
MEDIUM_NCELLS = 2500

# How much slower the medium fixture is allowed to run relative to the small
# fixture, expressed as a multiple of the plain cell-count ratio. Real Seurat
# steps (SCTransform, PCA, neighbor graphs, clustering) scale roughly linearly
# to mildly superlinearly with cell count, so a few-fold safety margin above
# the linear expectation comfortably covers normal variance while still
# catching a quadratic-or-worse regression (e.g. the do.call()/SCTransform
# hang this test suite failed to catch: see workflow/scripts/downsample_clusters.R).
SAFETY_MULTIPLIER = 4
CELL_RATIO = MEDIUM_NCELLS / SMALL_NCELLS
MAX_ALLOWED_RATIO = CELL_RATIO * SAFETY_MULTIPLIER

# This test uses nDownsampleReplicates=1 (unlike the other downsampling tests,
# which use 2) to keep the scaling comparison itself fast, so it can't reuse
# utils.assert_stability_table (which hardcodes EXPECTED_BOOTSTRAPS={1, 2}).
def _assert_basic_stability_table(path):
    columns, rows = read_tsv(path)
    assert rows, f"{path}: downsample output is empty"
    assert columns == EXPECTED_COLUMNS
    for row in rows:
        assert row["clusterid"] != ""
        assert int(row["bootstrap_number"]) == 1
        assert 0 <= float(row["max_jaccard"]) <= 1


def _run_downsample_clusters(sample, tmp_path, pytestconfig, root):
    seurat_dir = root / "testdata" / "downsampling" / "seurat_objects"
    results_dir = tmp_path / sample / "results"
    target = results_dir / f"{sample}_clusterdownsampling.tsv"

    cmd = [
        snakemake_executable(),
        str(target),
        "--snakefile",
        "workflow/Snakefile",
        "--configfile",
        "config/config.yaml",
        "--config",
        "workflow_mode=downsample_only",
        f"downsampleSeuratObjectDir={seurat_dir.as_posix()}",
        f"downsampleResultsDir={results_dir.as_posix()}",
        "nDownsampleReplicates=1",
        "downsampleRate=0.5",
        "workflowSeed=12345",
        "--profile",
        "none",
        "--workflow-profile",
        "none",
        "--executor",
        "local",
        "--cores",
        "1",
        "--jobs",
        "1",
        "--latency-wait",
        "30",
        "--rerun-incomplete",
        "--use-conda",
    ]
    conda_prefix = pytestconfig.getoption("--snakemake-conda-prefix")
    if conda_prefix:
        cmd.extend(["--conda-prefix", conda_prefix])

    t0 = time.monotonic()
    result = run_command(cmd, root, timeout=1800)
    elapsed = time.monotonic() - t0

    assert result.returncode == 0, combined_output(result)
    assert target.exists(), f"missing rule output: {target}"
    _assert_basic_stability_table(target)
    return elapsed


def test_downsample_clusters_runtime_scales_reasonably_with_ncells(tmp_path, pytestconfig):
    if not pytestconfig.getoption("--run-downsample-scaling"):
        pytest.skip("use --run-downsample-scaling to execute the downsample_clusters scaling check")

    root = repo_root()

    small_elapsed = _run_downsample_clusters(SMALL_SAMPLE, tmp_path, pytestconfig, root)
    medium_elapsed = _run_downsample_clusters(MEDIUM_SAMPLE, tmp_path, pytestconfig, root)

    ratio = medium_elapsed / small_elapsed
    assert ratio <= MAX_ALLOWED_RATIO, (
        f"downsample_clusters took {ratio:.1f}x longer on {MEDIUM_NCELLS} cells than on "
        f"{SMALL_NCELLS} cells ({small_elapsed:.1f}s -> {medium_elapsed:.1f}s), exceeding the "
        f"{MAX_ALLOWED_RATIO:.1f}x budget ({SAFETY_MULTIPLIER}x the {CELL_RATIO:.1f}x cell-count "
        "ratio). This suggests a non-linear performance regression in the R script."
    )
