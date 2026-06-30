import pytest

from utils import (
    assert_stability_table,
    combined_output,
    repo_root,
    run_command,
    snakemake_executable,
)


SAMPLE = "filtered_seurat_tenx_test"


def test_bootstrap_clusters_rule_produces_stability_table(tmp_path, pytestconfig):
    if not pytestconfig.getoption("--run-downsample-rule"):
        pytest.skip("use --run-downsample-rule to execute the downsample_clusters rule")

    root = repo_root()
    seurat_dir = root / "testdata" / "downsampling" / "seurat_objects"
    results_dir = tmp_path / "results"
    target = results_dir / f"{SAMPLE}_clusterdownsampling.tsv"

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
        "nDownsampleReplicates=2",
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

    result = run_command(cmd, root, timeout=1800)
    assert result.returncode == 0, combined_output(result)
    assert target.exists(), f"missing rule output: {target}"
    assert_stability_table(target)
