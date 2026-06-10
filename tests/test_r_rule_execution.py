import os
import shutil
import subprocess
from pathlib import Path


TEST_SAMPLE_SHEET = Path("testdata/samplesheet_test.tsv")
VALIDATOR = Path(__file__).with_name("validate_r_rule_outputs.R")


def run_command(cmd, repo_root, env=None, timeout=None):
    return subprocess.run(
        cmd,
        cwd=repo_root,
        env=env,
        text=True,
        capture_output=True,
        check=False,
        timeout=timeout,
    )


def combined_output(result):
    return result.stdout + result.stderr


def find_rscript_with_seuratobject(repo_root):
    candidates = []
    path_rscript = shutil.which("Rscript")
    if path_rscript:
        candidates.append(Path(path_rscript))
    candidates.extend(sorted((repo_root / ".snakemake" / "conda").glob("*/bin/Rscript")))

    seen = set()
    for candidate in candidates:
        if candidate in seen or not candidate.exists():
            continue
        seen.add(candidate)
        result = subprocess.run(
            [
                str(candidate),
                "-e",
                "quit(status = ifelse(requireNamespace('SeuratObject', quietly = TRUE), 0, 1))",
            ],
            text=True,
            capture_output=True,
            check=False,
        )
        if result.returncode == 0:
            return candidate

    raise AssertionError("could not find an Rscript with the SeuratObject package installed")


def test_tenx_marker_rule_chain_executes_r_scripts_and_produces_valid_outputs(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    sample_sheet = repo_root / TEST_SAMPLE_SHEET
    assert sample_sheet.exists(), f"missing test sample sheet: {sample_sheet}"
    assert VALIDATOR.exists(), f"missing R output validator: {VALIDATOR}"

    snakemake = shutil.which("snakemake")
    assert snakemake is not None, "snakemake is not available on PATH"

    results_dir = tmp_path / "results"
    seurat_rds = results_dir / "seurat_filtered" / "filtered_seurat_tenx_test.rds"
    marker_csv = results_dir / "seurat_filtered" / "filtered_seurat_tenx_test_markergenes.csv"

    env = os.environ.copy()
    env["SCRNASEQ_PREPROCESS_SEED"] = "12345"

    conda_frontend = os.environ.get("SNAKEMAKE_CONDA_FRONTEND")
    if conda_frontend is None and shutil.which("mamba"):
        conda_frontend = "mamba"

    cmd = [
        snakemake,
        str(marker_csv),
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
        "--use-conda",
        "--show-failed-logs",
        "--snakefile",
        "workflow/Snakefile",
        "--configfile",
        "config/config.yaml",
        "--config",
        f"sampleTable={TEST_SAMPLE_SHEET.as_posix()}",
        f"resultsDir={results_dir.as_posix()}",
    ]
    if conda_frontend:
        cmd.extend(["--conda-frontend", conda_frontend])

    result = run_command(cmd, repo_root, env=env, timeout=600)
    assert result.returncode == 0, combined_output(result)
    assert seurat_rds.exists(), f"missing Seurat RDS: {seurat_rds}"
    assert marker_csv.exists(), f"missing marker CSV: {marker_csv}"

    rscript = find_rscript_with_seuratobject(repo_root)
    validation = run_command(
        [str(rscript), str(VALIDATOR), str(seurat_rds), str(marker_csv)],
        repo_root,
        env=env,
        timeout=120,
    )
    assert validation.returncode == 0, combined_output(validation)
