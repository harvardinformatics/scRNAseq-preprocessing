import pytest

from utils import assert_stability_table, repo_root, run_downsample_clusters_rule


SAMPLE = "filtered_seurat_tenx_test"


def test_downsample_clusters_rule_produces_stability_table(tmp_path, pytestconfig):
    if not pytestconfig.getoption("--run-downsample-rule"):
        pytest.skip("use --run-downsample-rule to execute the downsample_clusters rule")

    root = repo_root()
    target = run_downsample_clusters_rule(SAMPLE, tmp_path, pytestconfig, root)
    assert_stability_table(target)
