import shutil
import subprocess
from pathlib import Path


TEST_SAMPLE_SHEET = Path("testdata/samplesheet_test.tsv")


def test_snakemake_lint_has_no_findings():
    repo_root = Path(__file__).resolve().parents[1]
    snakemake = shutil.which("snakemake")
    assert snakemake is not None, "snakemake is not available on PATH"

    cmd = [
        snakemake,
        "--lint",
        "--snakefile",
        "workflow/Snakefile",
        "--configfile",
        "config/config.yaml",
        "--config",
        f"sampleTable={TEST_SAMPLE_SHEET.as_posix()}",
        "resultsDir=testdata/results",
        "downsampleResultsDir=testdata/results/downsampling",
    ]
    result = subprocess.run(
        cmd,
        cwd=repo_root,
        text=True,
        capture_output=True,
        check=False,
        timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr
