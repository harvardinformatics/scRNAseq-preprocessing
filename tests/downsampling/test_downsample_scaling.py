import time

import pytest

from utils import assert_stability_table, repo_root, run_downsample_clusters_rule


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


def _timed_run(sample, tmp_path, pytestconfig, root):
    t0 = time.monotonic()
    target = run_downsample_clusters_rule(sample, tmp_path, pytestconfig, root, n_replicates=1)
    elapsed = time.monotonic() - t0
    assert_stability_table(target, expected_bootstraps={1})
    return elapsed


def test_downsample_clusters_runtime_scales_reasonably_with_ncells(tmp_path, pytestconfig):
    if not pytestconfig.getoption("--run-downsample-scaling"):
        pytest.skip("use --run-downsample-scaling to execute the downsample_clusters scaling check")

    root = repo_root()

    small_elapsed = _timed_run(SMALL_SAMPLE, tmp_path, pytestconfig, root)
    medium_elapsed = _timed_run(MEDIUM_SAMPLE, tmp_path, pytestconfig, root)

    ratio = medium_elapsed / small_elapsed
    assert ratio <= MAX_ALLOWED_RATIO, (
        f"downsample_clusters took {ratio:.1f}x longer on {MEDIUM_NCELLS} cells than on "
        f"{SMALL_NCELLS} cells ({small_elapsed:.1f}s -> {medium_elapsed:.1f}s), exceeding the "
        f"{MAX_ALLOWED_RATIO:.1f}x budget ({SAFETY_MULTIPLIER}x the {CELL_RATIO:.1f}x cell-count "
        "ratio). This suggests a non-linear performance regression in the R script."
    )
