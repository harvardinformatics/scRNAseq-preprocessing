import os
import shutil
import subprocess
from pathlib import Path


TEST_SAMPLE_SHEET = Path("testdata/samplesheet_test.tsv")


FAKE_CELLBENDER = r"""#!/usr/bin/env bash
set -euo pipefail

if [[ "${1:-}" == "remove-background" && "${2:-}" == "--help" ]]; then
    printf '%s\\n' "cellbender remove-background" "  --seed INTEGER"
    exit 0
fi

if [[ "${1:-}" != "remove-background" ]]; then
    echo "unexpected cellbender invocation: $*" >&2
    exit 2
fi
shift

input=""
output=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)
            input="$2"
            shift 2
            ;;
        --output)
            output="$2"
            shift 2
            ;;
        --seed)
            shift 2
            ;;
        --cuda)
            shift
            ;;
        *)
            shift
            ;;
    esac
done

if [[ -z "${input}" ]]; then
    echo "missing --input" >&2
    exit 3
fi
if [[ -z "${output}" ]]; then
    echo "missing --output" >&2
    exit 3
fi
if [[ ! -s "${input}" ]]; then
    echo "input does not exist or is empty: ${input}" >&2
    exit 4
fi

mkdir -p "$(dirname "${output}")"
printf 'fake raw cellbender output for %s\\n' "${input}" > "${output}"
filtered="${output%.h5}_filtered.h5"
printf 'fake filtered cellbender output for %s\\n' "${input}" > "${filtered}"
report="${output%.h5}_report.html"
printf '<html><body><h2>Automated assessment</h2><h2>Summary</h2><p>This learning curve looks normal.</p></body></html>\\n' > "${report}"
"""


def write_fake_cellbender(tmp_path):
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_cellbender = fake_bin / "cellbender"
    fake_cellbender.write_text(FAKE_CELLBENDER)
    fake_cellbender.chmod(0o755)
    return fake_bin


def test_cellbender_rule_uses_cellbender_filtered_output_convention(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    sample_sheet = repo_root / TEST_SAMPLE_SHEET
    assert sample_sheet.exists(), f"missing test sample sheet: {sample_sheet}"

    snakemake = shutil.which("snakemake")
    assert snakemake is not None, "snakemake is not available on PATH"

    fake_bin = write_fake_cellbender(tmp_path)
    results_dir = tmp_path / "results"
    base_output = results_dir / "cellbender" / "cellbender_test.h5"
    filtered_output = results_dir / "cellbender" / "cellbender_test_filtered.h5"
    report_output = results_dir / "cellbender" / "cellbender_test_report.html"
    status_output = results_dir / "cellbender" / "cellbender_test_adaptive_status.txt"

    env = os.environ.copy()
    env["PATH"] = f"{fake_bin}{os.pathsep}{env.get('PATH', '')}"

    cmd = [
        snakemake,
        str(base_output),
        str(filtered_output),
        "--profile",
        "none",
        "--workflow-profile",
        "none",
        "--executor",
        "local",
        "--cores",
        "1",
        "--latency-wait",
        "5",
        "--snakefile",
        "workflow/Snakefile",
        "--configfile",
        "config/config.yaml",
        "--config",
        f"sampleTable={TEST_SAMPLE_SHEET.as_posix()}",
        f"resultsDir={results_dir.as_posix()}",
    ]

    result = subprocess.run(
        cmd,
        cwd=repo_root,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert base_output.read_text().startswith("fake raw cellbender output")
    assert filtered_output.read_text().startswith("fake filtered cellbender output")
    # The adaptive wrapper also produces the report and per-sample status; a normal-looking
    # learning curve means the single initial run is kept, with no re-run.
    assert "This learning curve looks normal" in report_output.read_text()
    status = status_output.read_text()
    assert "outcome\tNORMAL_FIRST_TRY" in status
    assert "reran\tfalse" in status
