from utils import (
    TEST_RESULTS_DIR,
    TEST_SAMPLES,
    TEST_SEURAT_DIR,
    base_snakemake_cmd,
    combined_output,
    repo_root,
    run_command,
)


def test_testdata_fixtures_exist():
    root = repo_root()
    for sample in TEST_SAMPLES:
        assert (root / TEST_SEURAT_DIR / f"{sample}.rds").exists()


def test_testdata_dry_run_builds_expected_dag():
    root = repo_root()
    cmd = base_snakemake_cmd() + ["-np"]

    result = run_command(cmd, root)
    output = combined_output(result)

    assert result.returncode == 0, output
    assert "bootstrap_clusters" in output
    for sample in TEST_SAMPLES:
        expected = TEST_RESULTS_DIR / f"{sample}_clusterdownsampling.tsv"
        assert expected.as_posix() in output
