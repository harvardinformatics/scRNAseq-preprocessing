"""Parse-time preflight cell-count gate.

The `test` sample's CellRanger filtered matrix has 300 cells, so a threshold above 300
exercises the "below threshold" behavior and a threshold at/below 300 (the default 100)
passes it. All checks run as `snakemake -np` dry runs in workflow_mode=preprocess so that
skipping the only sample simply yields an empty DAG (in preprocess_and_downsample mode an
emptied sample set trips the separate "no inputs for downsampling" guard instead)."""

import shutil
import subprocess
from pathlib import Path


def run_dry_run(repo_root, results_dir, *config_overrides):
    snakemake = shutil.which("snakemake")
    assert snakemake is not None, "snakemake is not available on PATH"
    cmd = [
        snakemake,
        "-np",
        "--profile", "none",
        "--workflow-profile", "none",
        "--snakefile", "workflow/Snakefile",
        "--configfile", "config/config.yaml",
        "--config",
        "sampleTable=testdata/samplesheet_test.tsv",
        f"resultsDir={results_dir}",
        "workflow_mode=preprocess",
        *config_overrides,
    ]
    return subprocess.run(
        cmd, cwd=repo_root, text=True, capture_output=True, check=False, timeout=120
    )


def combined_output(result):
    return result.stdout + result.stderr


def test_preflight_passes_sample_above_threshold(tmp_path):
    # 300-cell test sample vs the default preflight_min_cells=100 -> runs normally.
    repo_root = Path(__file__).resolve().parents[1]
    result = run_dry_run(repo_root, tmp_path / "results")
    output = combined_output(result)

    assert result.returncode == 0, output
    assert "tenx2seuratrds" in output
    assert "below preflight_min_cells" not in output


def test_preflight_skip_drops_low_cell_sample(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    result = run_dry_run(
        repo_root, tmp_path / "results",
        "preflight_min_cells=500", "preflight_mode=skip",
    )
    output = combined_output(result)

    assert result.returncode == 0, output
    assert "below preflight_min_cells=500" in output
    assert "test (300 cells)" in output
    assert "skipping 1 low-cell sample(s)" in output
    assert "tenx2seuratrds" not in output  # the only sample was dropped from the DAG


def test_preflight_error_mode_hard_stops(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    result = run_dry_run(
        repo_root, tmp_path / "results",
        "preflight_min_cells=500", "preflight_mode=error",
    )
    output = combined_output(result)

    assert result.returncode != 0
    assert "fewer than preflight_min_cells=500" in output
    assert "test (300 cells)" in output


def test_preflight_warn_mode_reports_but_keeps_sample(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    result = run_dry_run(
        repo_root, tmp_path / "results",
        "preflight_min_cells=500", "preflight_mode=warn",
    )
    output = combined_output(result)

    assert result.returncode == 0, output
    assert "below preflight_min_cells=500" in output
    assert "test (300 cells)" in output
    assert "skipping" not in output
    assert "tenx2seuratrds" in output  # kept despite being below threshold


def test_preflight_off_mode_disables_check(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    result = run_dry_run(
        repo_root, tmp_path / "results",
        "preflight_min_cells=500", "preflight_mode=off",
    )
    output = combined_output(result)

    assert result.returncode == 0, output
    assert "[preflight]" not in output
    assert "tenx2seuratrds" in output
