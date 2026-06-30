import shutil

import pytest

from utils import (
    TEST_RESULTS_DIR,
    TEST_SAMPLES,
    assert_stability_table,
    assert_workflow_outputs_match_references,
    combined_output,
    expected_output_for_sample,
    repo_root,
    run_command,
    workflow_run_cmd,
)


def test_testdata_workflow_run(pytestconfig):
    if not pytestconfig.getoption("--run-downsample-workflow"):
        pytest.skip("use --run-downsample-workflow to execute the full downsampling workflow")

    root = repo_root()
    results_dir = root / TEST_RESULTS_DIR
    if results_dir.exists():
        shutil.rmtree(results_dir)

    cmd = workflow_run_cmd()
    conda_prefix = pytestconfig.getoption("--snakemake-conda-prefix")
    if conda_prefix:
        cmd.extend(["--conda-prefix", conda_prefix])

    result = run_command(cmd, root, timeout=3600)
    assert result.returncode == 0, combined_output(result)

    for sample in TEST_SAMPLES:
        output = root / expected_output_for_sample(sample)
        assert output.exists(), f"missing workflow output: {output}"
        assert_stability_table(output)

    assert_workflow_outputs_match_references(root)
