"""Check that the R environment Snakemake builds for the rules holds exactly the packages the pin file lists.

Snakemake creates the environment with `conda create --file workflow/envs/scrnaseq_r.linux-64.pin.txt`. If
that command fails, Snakemake logs a warning and builds the environment from `workflow/envs/scrnaseq_r.yml`
instead, which takes the newest build of every package that file leaves open. The reference outputs were
made with the pinned packages, so such an environment shows up only as small numeric differences in the full
workflow test, and pytest hides the warning because it captures Snakemake's output. This test fails on its
own when the environment differs from the pin file, and shows the error conda gave.
"""

import json
import os
import platform
import re
import shutil
import subprocess
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
ENV_FILE = Path("workflow/envs/scrnaseq_r.yml")
PIN_FILE = Path("workflow/envs/scrnaseq_r.linux-64.pin.txt")
TEST_SAMPLE_SHEET = Path("testdata/samplesheet_test.tsv")
FALLBACK_WARNING = "Failed to install conda environment from pin file"
SHOWN_DIFFERENCES = 100


def pinned_packages(pin_text):
    """File names of the packages an @EXPLICIT pin file lists (the end of each URL, without any #hash)."""
    packages = set()
    for line in pin_text.splitlines():
        line = line.strip()
        if line and not line.startswith(("#", "@")):
            packages.add(line.split("#", 1)[0].rsplit("/", 1)[-1])
    return packages


def installed_packages(env_dir):
    """File names of the packages conda records as installed in env_dir."""
    return {
        json.loads(record.read_text())["fn"]
        for record in (Path(env_dir) / "conda-meta").glob("*.json")
    }


def split_package_file(file_name):
    """(name, "version-build") of a package file name such as r-seurat-5.5.1-r45h1234_0.conda."""
    stem = re.sub(r"\.(conda|tar\.bz2)$", "", file_name)
    name, version, build = stem.rsplit("-", 2)
    return name, f"{version}-{build}"


def package_differences(pinned, installed):
    """One line per package that is installed in another build than the pinned one, or on one side only."""
    pinned_only = dict(split_package_file(name) for name in pinned - installed)
    installed_only = dict(split_package_file(name) for name in installed - pinned)
    return [
        f"{name}: pinned {pinned_only.get(name, '(none)')}, installed {installed_only.get(name, '(none)')}"
        for name in sorted(set(pinned_only) | set(installed_only))
    ]


def pin_install_failures(text):
    """Snakemake's warnings that installing from the pin file failed, each with the error conda gave."""
    lines = text.splitlines()
    failures = []
    for start, line in enumerate(lines):
        if FALLBACK_WARNING not in line:
            continue
        block = [line]
        for following in lines[start + 1 : start + 80]:
            # Snakemake logs this next, as it starts again from the environment file.
            if following.startswith("Downloading and installing remote packages"):
                break
            block.append(following)
        failures.append("\n".join(block))
    return failures


def snakemake_command(snakemake, results_dir, *options):
    cmd = [
        snakemake,
        *options,
        "--use-conda",
        "--profile",
        "none",
        "--workflow-profile",
        "none",
        "--snakefile",
        "workflow/Snakefile",
        "--configfile",
        "config/config.yaml",
    ]
    conda_frontend = os.environ.get("SNAKEMAKE_CONDA_FRONTEND")
    if conda_frontend is None and shutil.which("mamba"):
        conda_frontend = "mamba"
    if conda_frontend:
        cmd.extend(["--conda-frontend", conda_frontend])
    cmd.extend(
        [
            "--config",
            f"sampleTable={TEST_SAMPLE_SHEET.as_posix()}",
            f"resultsDir={results_dir.as_posix()}",
        ]
    )
    return cmd


def run_snakemake(cmd, timeout):
    result = subprocess.run(
        cmd, cwd=REPO_ROOT, text=True, capture_output=True, check=False, timeout=timeout
    )
    return result.returncode, result.stdout + result.stderr


def logged_pin_install_failures():
    """Pin-file warnings in the Snakemake logs under .snakemake/log, newest log first (at most two)."""
    logs = sorted(
        (REPO_ROOT / ".snakemake" / "log").glob("*.snakemake.log"),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    found = []
    for log in logs:
        for failure in pin_install_failures(log.read_text(errors="replace")):
            found.append(f"(from {log.relative_to(REPO_ROOT)})\n{failure}")
        if len(found) >= 2:
            break
    return found[:2]


def test_package_differences_name_each_package_that_differs(tmp_path):
    pin_text = "\n".join(
        [
            "# platform: linux-64",
            "@EXPLICIT",
            "https://conda.anaconda.org/conda-forge/linux-64/r-base-4.5.3-h1_0.conda",
            "https://conda.anaconda.org/conda-forge/noarch/r-seurat-5.5.1-r45h2_0.conda#0123abcd",
            "https://conda.anaconda.org/conda-forge/linux-64/libfoo-1.0-h3_0.tar.bz2",
        ]
    )
    conda_meta = tmp_path / "conda-meta"
    conda_meta.mkdir()
    for file_name in ["r-base-4.5.3-h1_0.conda", "r-seurat-5.6.0-r45h9_0.conda", "libbar-2.0-h4_0.conda"]:
        record = conda_meta / (file_name.rsplit(".", 1)[0] + ".json")
        record.write_text(json.dumps({"fn": file_name}))

    assert package_differences(pinned_packages(pin_text), installed_packages(tmp_path)) == [
        "libbar: pinned (none), installed 2.0-h4_0",
        "libfoo: pinned 1.0-h3_0, installed (none)",
        "r-seurat: pinned 5.5.1-r45h2_0, installed 5.6.0-r45h9_0",
    ]


def test_pin_install_failures_keep_the_error_conda_gave():
    log = "\n".join(
        [
            "Creating conda environment workflow/envs/scrnaseq_r.yml...",
            "Using pinnings from workflow/envs/scrnaseq_r.linux-64.pin.txt.",
            "Downloading and installing remote packages.",
            f"{FALLBACK_WARNING} (workflow/envs/scrnaseq_r.linux-64.pin.txt). Trying regular environment definition file.",
            "Error message:",
            "CondaHTTPError: HTTP 404 NOT FOUND",
            "Downloading and installing remote packages.",
            "Environment for workflow/envs/scrnaseq_r.yml created (location: .snakemake/conda/0123_)",
        ]
    )

    assert pin_install_failures(log) == [
        f"{FALLBACK_WARNING} (workflow/envs/scrnaseq_r.linux-64.pin.txt). Trying regular environment definition file.\n"
        "Error message:\n"
        "CondaHTTPError: HTTP 404 NOT FOUND"
    ]
    assert pin_install_failures("Using pinnings from workflow/envs/scrnaseq_r.linux-64.pin.txt.") == []


def test_rule_r_environment_matches_the_pin_file(tmp_path, github_error):
    if platform.system() != "Linux" or platform.machine() != "x86_64":
        pytest.skip("the R environment is pinned for linux-64 only")

    snakemake = shutil.which("snakemake")
    assert snakemake is not None, "snakemake is not available on PATH"
    results_dir = tmp_path / "results"

    # Build the environment as the rules would (nothing to do when it exists), then ask Snakemake where it is.
    returncode, created = run_snakemake(
        snakemake_command(snakemake, results_dir, "--conda-create-envs-only", "--cores", "1"),
        timeout=1800,
    )
    assert returncode == 0, created
    returncode, listed = run_snakemake(
        snakemake_command(snakemake, results_dir, "--list-conda-envs"), timeout=600
    )
    assert returncode == 0, listed
    locations = [
        line.split("\t")[-1].strip()
        for line in listed.splitlines()
        if line.split("\t")[0].endswith(ENV_FILE.name)
    ]
    assert len(locations) == 1, f"expected one environment for {ENV_FILE}:\n{listed}"
    env_dir = REPO_ROOT / locations[0]
    assert (env_dir / "conda-meta").is_dir(), f"no conda environment at {env_dir}"

    differences = package_differences(
        pinned_packages((REPO_ROOT / PIN_FILE).read_text()), installed_packages(env_dir)
    )
    if not differences:
        return

    report = [f"The R environment {locations[0]} differs from {PIN_FILE} in {len(differences)} package(s):"]
    report += differences[:SHOWN_DIFFERENCES]
    if len(differences) > SHOWN_DIFFERENCES:
        report.append(f"... and {len(differences) - SHOWN_DIFFERENCES} more")
    # Snakemake copies the file it builds an environment from next to the environment.
    if Path(f"{env_dir}.yaml").exists():
        report.append(
            f"Snakemake built it from {ENV_FILE}, as it does after installing from the pin file fails."
        )
    failures = pin_install_failures(created) or logged_pin_install_failures()
    if failures:
        report.append("Snakemake's warning:")
        report += failures
    else:
        report.append(
            "No pin-file warning in this run's Snakemake output or in .snakemake/log; the environment may "
            "have been built earlier (for example, restored from the GitHub Actions cache of .snakemake/conda)."
        )

    github_error("R environment differs from the pin file", "\n".join(report))
    pytest.fail("\n".join(report), pytrace=False)
