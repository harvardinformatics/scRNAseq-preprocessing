"""Parse-time preflight cell-count gate.

The `test` sample's CellRanger filtered matrix has 300 cells, so a threshold above 300
exercises the "below threshold" behavior and a threshold at/below 300 (the default 100)
passes it. All checks run as `snakemake -np` dry runs in workflow_mode=preprocess so that
skipping the only sample simply yields an empty DAG (in preprocess_and_downsample mode an
emptied sample set trips the separate "no inputs for downsampling" guard instead)."""

import gzip
import shutil
import subprocess
from pathlib import Path


def run_dry_run(repo_root, results_dir, *config_overrides, sample_table="testdata/samplesheet_test.tsv"):
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
        f"sampleTable={sample_table}",
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


# --- CellBender raw-matrix preflight (integrity -> error, droplet sufficiency -> warn) -------

def _write_gz(path, text):
    with gzip.open(path, "wt") as handle:
        handle.write(text)


def make_datadir(base, filtered_cells, raw_features, raw_barcodes, raw_nnz):
    """A minimal 10x datadir: `filtered_cells` filtered barcodes and a raw matrix.mtx.gz whose
    header advertises the given (features, barcodes, nnz). Enough to pass validate_sample_sheet
    and exercise the CellBender raw-matrix checks. filtered_cells >= 100 so the cell-count gate
    doesn't skip the sample before the CellBender check runs."""
    filtered = base / "filtered_feature_bc_matrix"
    filtered.mkdir(parents=True)
    _write_gz(filtered / "barcodes.tsv.gz", "".join(f"CELL{i}-1\n" for i in range(filtered_cells)))

    raw = base / "raw_feature_bc_matrix"
    raw.mkdir(parents=True)
    body = "1 1 1\n" if raw_nnz > 0 else ""
    _write_gz(
        raw / "matrix.mtx.gz",
        "%%MatrixMarket matrix coordinate integer general\n"
        f"{raw_features} {raw_barcodes} {raw_nnz}\n" + body,
    )
    (base / "raw_feature_bc_matrix.h5").write_text("")  # existence only; not read by the checks
    return base


def write_sample_sheet(path, sampleid, datadir):
    path.write_text(f"sampleid\ttenx_datadir\n{sampleid}\t{datadir}\n")
    return path


def test_cellbender_preflight_warns_on_too_few_raw_droplets(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    datadir = make_datadir(tmp_path / "lowraw", filtered_cells=300,
                           raw_features=2000, raw_barcodes=400, raw_nnz=1000)  # 400 < 2*300
    sheet = write_sample_sheet(tmp_path / "sheet.tsv", "lowraw", datadir)
    result = run_dry_run(repo_root, tmp_path / "results", sample_table=sheet)
    output = combined_output(result)

    assert result.returncode == 0, output  # a warning is non-blocking
    assert "[preflight:cellbender] WARNING" in output
    assert "400 raw droplets vs 300 called cells" in output
    assert "tenx2seuratrds" in output  # the sample still runs


def test_cellbender_preflight_errors_on_empty_raw_matrix(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    datadir = make_datadir(tmp_path / "emptyraw", filtered_cells=300,
                           raw_features=0, raw_barcodes=0, raw_nnz=0)
    sheet = write_sample_sheet(tmp_path / "sheet.tsv", "emptyraw", datadir)
    result = run_dry_run(repo_root, tmp_path / "results", sample_table=sheet)
    output = combined_output(result)

    assert result.returncode != 0
    assert "[preflight:cellbender]" in output
    assert "empty/degenerate" in output


def test_cellbender_preflight_errors_when_raw_smaller_than_filtered(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    datadir = make_datadir(tmp_path / "mismatch", filtered_cells=300,
                           raw_features=2000, raw_barcodes=200, raw_nnz=1000)  # raw < filtered
    sheet = write_sample_sheet(tmp_path / "sheet.tsv", "mismatch", datadir)
    result = run_dry_run(repo_root, tmp_path / "results", sample_table=sheet)
    output = combined_output(result)

    assert result.returncode != 0
    assert "fewer barcodes (200) than the filtered matrix (300)" in output


def test_cellbender_preflight_passes_a_true_raw_matrix(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    datadir = make_datadir(tmp_path / "goodraw", filtered_cells=300,
                           raw_features=2000, raw_barcodes=20000, raw_nnz=5000)
    sheet = write_sample_sheet(tmp_path / "sheet.tsv", "goodraw", datadir)
    result = run_dry_run(repo_root, tmp_path / "results", sample_table=sheet)
    output = combined_output(result)

    assert result.returncode == 0, output
    assert "[preflight:cellbender]" not in output
    assert "tenx2seuratrds" in output
