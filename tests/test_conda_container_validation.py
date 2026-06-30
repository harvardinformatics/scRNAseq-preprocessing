import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml


ENV_DIR = Path("workflow/envs")
RULE_DIR = Path("workflow/rules")
TEST_SAMPLE_SHEET = Path("testdata/samplesheet_test.tsv")
EXPECTED_ENV_FILES = {
    "cellbender.yml",
    "doubletfinder.yml",
    "downsample_clusters.yml",
    "emptydrops.yml",
    "posthocfilter.yml",
    "scdblfinder.yml",
    "soupx.yml",
    "tenx2seuratrds.yml",
}
CONDA_REFERENCE_RE = re.compile(r"conda:\s*\n\s*['\"]\.\./envs/([^'\"]+)['\"]")
CONTAINER_RE = re.compile(r"container:\s*\n\s*['\"]([^'\"]+)['\"]")

R_IMPORTS_BY_ENV = {
    "cellbender.yml": ["Seurat", "tidyverse", "bluster"],
    "doubletfinder.yml": ["Seurat", "tidyverse", "remotes", "fields", "Matrix", "KernSmooth", "ROCR", "igraph", "glmGamPoi", "bluster"],
    "downsample_clusters.yml": ["Seurat", "tidyverse", "glmGamPoi"],
    "emptydrops.yml": ["Seurat", "tidyverse", "DropletUtils", "scater", "glmGamPoi", "bluster", "Matrix", "R.utils"],
    "posthocfilter.yml": ["Seurat", "tidyverse", "glmGamPoi", "scater", "bluster"],
    "scdblfinder.yml": ["Seurat", "tidyverse", "scDblFinder", "glmGamPoi", "bluster"],
    "soupx.yml": ["Seurat", "tidyverse", "glmGamPoi", "bluster", "SoupX"],
    "tenx2seuratrds.yml": ["Seurat", "tidyverse", "scCustomize", "hdf5r", "glmGamPoi", "bluster", "presto"],
}
PYTHON_IMPORTS_BY_ENV = {
    "cellbender.yml": ["cellbender"],
}


def repo_root():
    return Path(__file__).resolve().parents[1]


def read_yaml(path):
    with path.open() as handle:
        return yaml.safe_load(handle)


def combined_output(result):
    return result.stdout + result.stderr


def run_command(cmd, cwd, env=None, timeout=120):
    return subprocess.run(
        cmd,
        cwd=cwd,
        env=env,
        text=True,
        capture_output=True,
        check=False,
        timeout=timeout,
    )


def available_conda_frontend():
    requested = os.environ.get("SNAKEMAKE_CONDA_FRONTEND")
    if requested and shutil.which(requested):
        return shutil.which(requested)
    for candidate in ("mamba", "conda"):
        path = shutil.which(candidate)
        if path:
            return path
    raise AssertionError("neither mamba nor conda is available on PATH")


def r_require_namespace_expr(packages):
    package_vector = ", ".join(repr(package) for package in packages)
    return (
        f"packages <- c({package_vector}); "
        "missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]; "
        "if (length(missing)) { stop(paste('missing R packages:', paste(missing, collapse = ',')), call. = FALSE) }"
    )


def test_selected_conda_env_name_is_known(pytestconfig):
    selected = pytestconfig.getoption("--conda-env-name")
    if selected is not None:
        assert selected in EXPECTED_ENV_FILES, f"unknown workflow env requested: {selected}"


def test_workflow_conda_env_files_are_well_formed_and_referenced_envs_exist():
    root = repo_root()
    env_paths = sorted((root / ENV_DIR).glob("*.yml"))
    observed = {path.name for path in env_paths}
    assert observed == EXPECTED_ENV_FILES

    for env_path in env_paths:
        env = read_yaml(env_path)
        assert env["channels"] == ["conda-forge", "bioconda"], env_path
        assert env.get("channel_priority") == "strict", env_path
        assert env.get("dependencies"), env_path

    referenced = set()
    for rule_path in sorted((root / RULE_DIR).glob("*.smk")):
        referenced.update(CONDA_REFERENCE_RE.findall(rule_path.read_text()))

    missing_references = sorted(name for name in referenced if not (root / ENV_DIR / name).exists())
    assert not missing_references, "rule conda env references are missing: " + ", ".join(missing_references)


def test_workflow_container_declarations_are_explicit_and_recognized():
    root = repo_root()
    containers = []
    for rule_path in sorted((root / RULE_DIR).glob("*.smk")):
        for uri in CONTAINER_RE.findall(rule_path.read_text()):
            containers.append((rule_path, uri))

    assert containers, "no workflow container declarations found"
    assert containers == [(root / "workflow/rules/cellbender.smk", "docker://us.gcr.io/broad-dsde-methods/cellbender:latest")]
    for _, uri in containers:
        assert uri.startswith("docker://")
        assert ":" in uri.removeprefix("docker://"), f"container URI is missing a tag: {uri}"


@pytest.mark.parametrize("env_name", sorted(EXPECTED_ENV_FILES))
def test_conda_env_solves_and_key_packages_import(tmp_path, pytestconfig, env_name):
    if not pytestconfig.getoption("--run-conda-validation"):
        pytest.skip("use --run-conda-validation to create and import-check workflow conda envs")
    selected = pytestconfig.getoption("--conda-env-name")
    if selected is not None and env_name != selected:
        pytest.skip(f"only validating requested workflow env: {selected}")

    root = repo_root()
    env_path = root / ENV_DIR / env_name
    prefix = tmp_path / env_name.removesuffix(".yml")
    conda = available_conda_frontend()

    create = run_command(
        [conda, "env", "create", "--yes", "--prefix", str(prefix), "--file", str(env_path)],
        root,
        timeout=1800,
    )
    assert create.returncode == 0, combined_output(create)

    r_packages = R_IMPORTS_BY_ENV.get(env_name, [])
    if r_packages:
        rscript = prefix / "bin" / "Rscript"
        assert rscript.exists(), f"missing Rscript in {prefix}"
        imports = run_command([str(rscript), "-e", r_require_namespace_expr(r_packages)], root, timeout=300)
        assert imports.returncode == 0, combined_output(imports)

    python_packages = PYTHON_IMPORTS_BY_ENV.get(env_name, [])
    if python_packages:
        python = prefix / "bin" / "python"
        assert python.exists(), f"missing python in {prefix}"
        statements = "; ".join(f"import {package}" for package in python_packages)
        imports = run_command([str(python), "-c", statements], root, timeout=300)
        assert imports.returncode == 0, combined_output(imports)


def test_doubletfinder_github_action_exports_token_for_remotes():
    root = repo_root()
    workflow = read_yaml(root / ".github/workflows/tests.yml")
    job = workflow["jobs"]["doubletfinder-install"]

    assert job.get("permissions", {}).get("contents") == "read"
    assert job.get("env", {}).get("GITHUB_PAT") == "${{ github.token }}"


def test_cellbender_container_can_be_pulled(tmp_path, pytestconfig):
    if not pytestconfig.getoption("--run-container-validation"):
        pytest.skip("use --run-container-validation to pull workflow containers")

    root = repo_root()
    uri = "docker://us.gcr.io/broad-dsde-methods/cellbender:latest"
    docker_image = uri.removeprefix("docker://")

    if shutil.which("docker"):
        result = run_command(["docker", "pull", docker_image], root, timeout=1800)
    elif shutil.which("apptainer"):
        result = run_command(["apptainer", "pull", str(tmp_path / "cellbender.sif"), uri], root, timeout=1800)
    elif shutil.which("singularity"):
        result = run_command(["singularity", "pull", str(tmp_path / "cellbender.sif"), uri], root, timeout=1800)
    else:
        pytest.skip("docker, apptainer, or singularity is required for container validation")

    assert result.returncode == 0, combined_output(result)


def test_doubletfinder_github_install_rule_executes(tmp_path, pytestconfig):
    if not pytestconfig.getoption("--run-doubletfinder-install"):
        pytest.skip("use --run-doubletfinder-install to test the networked GitHub install rule")

    root = repo_root()
    sample_sheet = root / TEST_SAMPLE_SHEET
    assert sample_sheet.exists(), f"missing test sample sheet: {sample_sheet}"

    snakemake = shutil.which("snakemake")
    assert snakemake is not None, "snakemake is not available on PATH"

    results_dir = tmp_path / "results"
    conda_prefix = tmp_path / "snakemake-conda"
    target = results_dir / "doubletfinder_installed.txt"
    env = os.environ.copy()
    conda_frontend = os.environ.get("SNAKEMAKE_CONDA_FRONTEND")
    if conda_frontend is None and shutil.which("mamba"):
        conda_frontend = "mamba"

    cmd = [
        snakemake,
        str(target),
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
        "--conda-prefix",
        str(conda_prefix),
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

    result = run_command(cmd, root, env=env, timeout=1800)
    assert result.returncode == 0, combined_output(result)
    assert target.exists(), f"missing DoubletFinder install sentinel: {target}"

    rscript_candidates = sorted(conda_prefix.glob("*/bin/Rscript"))
    assert rscript_candidates, f"no Snakemake conda Rscript found under {conda_prefix}"
    import_result = None
    for rscript in rscript_candidates:
        import_result = run_command(
            [str(rscript), "-e", "quit(status = ifelse(requireNamespace('DoubletFinder', quietly = TRUE), 0, 1))"],
            root,
            timeout=120,
        )
        if import_result.returncode == 0:
            break
    assert import_result is not None and import_result.returncode == 0, combined_output(import_result)
