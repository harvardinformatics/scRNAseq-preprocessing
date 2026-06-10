import shutil
import subprocess
from pathlib import Path


REQUIRED_TENX_CHILDREN = [
    "filtered_feature_bc_matrix",
    "raw_feature_bc_matrix",
    "raw_feature_bc_matrix.h5",
]


def run_snakemake_with_sample_sheet(repo_root, sample_sheet, results_dir, *extra_args):
    snakemake = shutil.which("snakemake")
    assert snakemake is not None, "snakemake is not available on PATH"

    cmd = [
        snakemake,
        *extra_args,
        "--profile",
        "none",
        "--workflow-profile",
        "none",
        "--snakefile",
        "workflow/Snakefile",
        "--configfile",
        "config/config.yaml",
        "--config",
        f"sampleTable={sample_sheet}",
        f"resultsDir={results_dir}",
    ]
    return subprocess.run(
        cmd,
        cwd=repo_root,
        text=True,
        capture_output=True,
        check=False,
        timeout=120,
    )


def combined_output(result):
    return result.stdout + result.stderr


def write_sample_sheet(path, rows, header="sampleid\ttenx_datadir"):
    lines = [header]
    lines.extend(f"{sampleid}\t{tenx_datadir}" for sampleid, tenx_datadir in rows)
    path.write_text("\n".join(lines) + "\n")
    return path


def test_malformed_sample_sheet_missing_required_column_fails(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    sample_sheet = tmp_path / "malformed_samplesheet.tsv"
    sample_sheet.write_text("sampleid\nmissing_tenx_column\n")

    result = run_snakemake_with_sample_sheet(
        repo_root,
        sample_sheet,
        tmp_path / "results",
        "-np",
    )
    output = combined_output(result)

    assert result.returncode != 0
    assert "Invalid sample sheet" in output
    assert "missing required column(s): tenx_datadir" in output


def test_duplicate_sample_ids_fail_validation(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    sample_sheet = write_sample_sheet(
        tmp_path / "duplicate_samplesheet.tsv",
        [("duplicate", "testdata"), ("duplicate", "testdata")],
    )

    result = run_snakemake_with_sample_sheet(
        repo_root,
        sample_sheet,
        tmp_path / "results",
        "-np",
    )
    output = combined_output(result)

    assert result.returncode != 0
    assert "Invalid sample sheet" in output
    assert "duplicate sampleid value(s): duplicate" in output


def test_missing_tenx_datadir_fails_validation(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    missing_data_dir = tmp_path / "does_not_exist"
    sample_sheet = write_sample_sheet(
        tmp_path / "missing_data_samplesheet.tsv",
        [("missing_path", missing_data_dir)],
    )

    result = run_snakemake_with_sample_sheet(
        repo_root,
        sample_sheet,
        tmp_path / "results",
        "-np",
    )
    output = combined_output(result)

    assert result.returncode != 0
    assert "Invalid sample sheet" in output
    assert "tenx_datadir does not exist" in output
    assert str(missing_data_dir) in output


def test_incomplete_tenx_datadir_fails_with_missing_required_children(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    incomplete_data_dir = tmp_path / "incomplete_testdata"
    incomplete_data_dir.mkdir()
    (incomplete_data_dir / "filtered_feature_bc_matrix").mkdir()

    sample_sheet = write_sample_sheet(
        tmp_path / "incomplete_data_samplesheet.tsv",
        [("incomplete", incomplete_data_dir)],
    )

    result = run_snakemake_with_sample_sheet(
        repo_root,
        sample_sheet,
        tmp_path / "results",
        "-np",
    )
    output = combined_output(result)

    assert result.returncode != 0
    assert "Invalid sample sheet" in output
    assert "missing raw_feature_bc_matrix" in output
    assert "missing raw_feature_bc_matrix.h5" in output


def test_absolute_tenx_datadir_is_normalized_and_used_for_rule_inputs(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    absolute_testdata = (repo_root / "testdata").resolve()
    sample_sheet = write_sample_sheet(
        tmp_path / "absolute_data_samplesheet.tsv",
        [("absolute_path", absolute_testdata)],
    )

    result = run_snakemake_with_sample_sheet(
        repo_root,
        sample_sheet,
        tmp_path / "results",
        "-np",
    )
    output = combined_output(result)

    assert result.returncode == 0, output
    for child in REQUIRED_TENX_CHILDREN:
        assert str(absolute_testdata / child) in output


def test_sample_sheet_with_multiple_samples_expands_outputs_for_each_sample(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    sample_sheet = write_sample_sheet(
        tmp_path / "multi_sample_samplesheet.tsv",
        [("test_a", "testdata"), ("test_b", "testdata")],
    )
    results_dir = tmp_path / "results"

    result = run_snakemake_with_sample_sheet(
        repo_root,
        sample_sheet,
        results_dir,
        "-np",
    )
    output = combined_output(result)

    assert result.returncode == 0, output
    for sample_id in ["test_a", "test_b"]:
        assert f"wildcards: sample={sample_id}" in output
        assert str(results_dir / "seurat_filtered" / f"filtered_seurat_tenx_{sample_id}_markergenes.csv") in output
        assert str(results_dir / "cellbender" / f"cellbender_{sample_id}.h5") in output
        assert str(results_dir / "cellbender" / f"cellbender_{sample_id}_filtered.h5") in output
