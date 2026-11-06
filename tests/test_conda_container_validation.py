import os
import platform
import re
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml


ENV_DIR = Path("workflow/envs")
RULE_DIR = Path("workflow/rules")
TEST_SAMPLE_SHEET = Path("testdata/samplesheet_test.tsv")
# Every R step runs in this one environment, so all methods share R and package versions.
SHARED_R_ENV = "scrnaseq_r.yml"
EXPECTED_ENV_FILES = {SHARED_R_ENV}
CONDA_REFERENCE_RE = re.compile(r"conda:\s*\n\s*['\"]\.\./envs/([^'\"]+)['\"]")
CONTAINER_RE = re.compile(r"container:\s*\n\s*['\"]([^'\"]+)['\"]")
CELLBENDER_CONTAINER_URI = "docker://us.gcr.io/broad-dsde-methods/cellbender@sha256:093f2caf1ce4acae4541ea45e52ab7b220ca131ec73b4d1f664b85fe12850bae"
DOUBLETFINDER_PINNED_RE = re.compile(r"chris-mcginnis-ucsf/DoubletFinder@[0-9a-f]{40}\b")

R_IMPORTS_BY_ENV = {
    SHARED_R_ENV: [
        "Seurat", "SeuratObject", "sctransform", "tidyverse", "glmGamPoi", "bluster", "Matrix", "R.utils",
        "RcppParallel", "DropletUtils", "rhdf5", "scater", "scDblFinder", "BiocParallel", "decontX", "SoupX",
        "scCustomize", "hdf5r", "presto", "remotes", "fields", "KernSmooth", "ROCR", "igraph",
    ],
}
PYTHON_IMPORTS_BY_ENV = {}


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


def linux_conda_subdir():
    if platform.system() == "Linux" and platform.machine() in {"x86_64", "AMD64"}:
        return "linux-64"
    return None


def pin_file_for_env(env_path):
    subdir = linux_conda_subdir()
    if subdir is None:
        return None
    return env_path.with_suffix(f".{subdir}.pin.txt")


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

        pin_path = pin_file_for_env(env_path)
        if pin_path is not None:
            assert pin_path.exists(), f"missing Snakemake conda pin file: {pin_path}"
            pin_text = pin_path.read_text()
            assert "@EXPLICIT" in pin_text, f"pin file is not an explicit conda spec: {pin_path}"
            assert "conda-forge" in pin_text or "bioconda" in pin_text, pin_path

    referenced = set()
    for rule_path in sorted((root / RULE_DIR).glob("*.smk")):
        referenced.update(CONDA_REFERENCE_RE.findall(rule_path.read_text()))

    missing_references = sorted(name for name in referenced if not (root / ENV_DIR / name).exists())
    assert not missing_references, "rule conda env references are missing: " + ", ".join(missing_references)


def test_every_conda_rule_uses_the_shared_r_environment():
    root = repo_root()
    rules_by_env = {}
    for rule_path in sorted((root / RULE_DIR).glob("*.smk")):
        for name in CONDA_REFERENCE_RE.findall(rule_path.read_text()):
            rules_by_env.setdefault(name, set()).add(rule_path.name)
    assert set(rules_by_env) == {SHARED_R_ENV}, rules_by_env


def test_doubletfinder_install_is_pinned_to_a_commit():
    rule_text = (repo_root() / RULE_DIR / "install_doubletfinder.smk").read_text()
    assert DOUBLETFINDER_PINNED_RE.search(rule_text), "install_doubletfinder must install a fixed DoubletFinder commit"


def test_workflow_container_declarations_are_explicit_and_recognized():
    root = repo_root()
    containers = []
    for rule_path in sorted((root / RULE_DIR).glob("*.smk")):
        for uri in CONTAINER_RE.findall(rule_path.read_text()):
            containers.append((rule_path, uri))

    assert containers, "no workflow container declarations found"
    # The `cellbender` rule and the `cellbender_adaptive_summary` rule (which reuses the same
    # image for a stdlib aggregation) are the only container declarations, and both pin the
    # identical CellBender digest.
    cellbender_smk = root / "workflow/rules/cellbender.smk"
    assert containers == [
        (cellbender_smk, CELLBENDER_CONTAINER_URI),
        (cellbender_smk, CELLBENDER_CONTAINER_URI),
    ]
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

    pin_path = pin_file_for_env(env_path)
    if pin_path is not None and pin_path.exists():
        create_cmd = [conda, "create", "--yes", "--prefix", str(prefix), "--file", str(pin_path)]
    else:
        create_cmd = [conda, "env", "create", "--yes", "--prefix", str(prefix), "--file", str(env_path)]

    create = run_command(
        create_cmd,
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
    uri = CELLBENDER_CONTAINER_URI
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
