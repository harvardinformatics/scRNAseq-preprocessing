"""Adaptive CellBender re-run wrapper (workflow/scripts/cellbender_adaptive_run.py) and the
run-level summary (cellbender_adaptive_summary.py).

`assess_report` is unit-tested directly. The orchestration is tested by driving `main()` against a
fake `cellbender` on PATH whose report verdict is scripted per run (initial vs re-run) through
environment variables, so every branch -- keep initial, re-run and resolve, re-run and fail to
resolve, adaptive off, and the hard error on an unrecognized report -- is exercised without a GPU."""

import importlib.util
import os
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"


def _load(module_name):
    path = SCRIPTS / f"{module_name}.py"
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


run_mod = _load("cellbender_adaptive_run")
summary_mod = _load("cellbender_adaptive_summary")


# A fake `cellbender`: advertises --seed, and on each real run writes the .h5/_filtered.h5 outputs
# plus a _report.html whose verdict is chosen by a per-invocation counter -- the first real run
# reads $FAKE_INITIAL_PHRASE, the second reads $FAKE_RERUN_PHRASE (values: normal|rerun|unknown).
FAKE_CELLBENDER = r"""#!/usr/bin/env bash
set -euo pipefail

if [[ "${1:-}" == "remove-background" && "${2:-}" == "--help" ]]; then
    printf '%s\n' "cellbender remove-background" "  --seed INTEGER"
    exit 0
fi
if [[ "${1:-}" != "remove-background" ]]; then
    echo "unexpected cellbender invocation: $*" >&2
    exit 2
fi
shift

output=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --output) output="$2"; shift 2 ;;
        --input|--learning-rate|--seed) shift 2 ;;
        *) shift ;;
    esac
done

n=1
if [[ -f "${FAKE_COUNTER}" ]]; then n=$(( $(cat "${FAKE_COUNTER}") + 1 )); fi
echo "${n}" > "${FAKE_COUNTER}"
if [[ "${n}" == "1" ]]; then phrase="${FAKE_INITIAL_PHRASE:-normal}"; else phrase="${FAKE_RERUN_PHRASE:-normal}"; fi

case "${phrase}" in
    normal)  summary="This learning curve looks normal." ;;
    rerun)   summary="Consider re-running with half the current learning rate to compare the results." ;;
    *)       summary="The assessment is inconclusive here." ;;
esac

mkdir -p "$(dirname "${output}")"
printf 'fake raw for run %s\n' "${n}" > "${output}"
printf 'fake filtered for run %s\n' "${n}" > "${output%.h5}_filtered.h5"
printf '<html><body><h2>Automated assessment</h2><h2>Summary</h2><p>%s</p></body></html>\n' \
    "${summary}" > "${output%.h5}_report.html"
"""


def _install_fake(tmp_path, monkeypatch, initial_phrase, rerun_phrase="normal"):
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake = fake_bin / "cellbender"
    fake.write_text(FAKE_CELLBENDER)
    fake.chmod(0o755)
    # Prepend so the fake wins as `cellbender` while bash/env/coreutils still resolve.
    monkeypatch.setenv("PATH", f"{fake_bin}:{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_COUNTER", str(tmp_path / "counter"))
    monkeypatch.setenv("FAKE_INITIAL_PHRASE", initial_phrase)
    monkeypatch.setenv("FAKE_RERUN_PHRASE", rerun_phrase)


def _run(tmp_path, adaptive="true"):
    """Drive main() and return (status_dict, results_dir Path)."""
    input_h5 = tmp_path / "raw.h5"
    input_h5.write_text("raw input")
    results = tmp_path / "results" / "cellbender"
    run_mod.main([
        "--input", str(input_h5),
        "--output-base", str(results / "cellbender_test.h5"),
        "--output-filtered", str(results / "cellbender_test_filtered.h5"),
        "--report", str(results / "cellbender_test_report.html"),
        "--status", str(results / "cellbender_test_adaptive_status.txt"),
        "--sample", "test",
        "--learning-rate", "0.0001",
        "--adaptive", adaptive,
        "--seed", "12345",
        "--workdir", str(tmp_path / "scratch"),
    ])
    status = {}
    for line in (results / "cellbender_test_adaptive_status.txt").read_text().splitlines():
        key, value = line.split("\t", 1)
        status[key] = value
    return status, results


# --- assess_report unit tests --------------------------------------------------------------------

def test_assess_report_normal_ignores_tags_and_case():
    html = "<div><h2>Summary</h2><p>This <b>learning</b> curve <i>looks</i> NORMAL.</p></div>"
    assert run_mod.assess_report(html) == "normal"


def test_assess_report_detects_rerun_suggestion():
    html = "<p>Consider re-running with half the current learning rate to compare.</p>"
    assert run_mod.assess_report(html) == "rerun"


def test_assess_report_unknown_when_neither_phrase_present():
    assert run_mod.assess_report("<p>Training finished. Elbo converged.</p>") == "unknown"


# --- orchestration -------------------------------------------------------------------------------

def test_normal_first_try_keeps_initial(tmp_path, monkeypatch):
    _install_fake(tmp_path, monkeypatch, initial_phrase="normal")
    status, results = _run(tmp_path)

    assert status["outcome"] == "NORMAL_FIRST_TRY"
    assert status["reran"] == "false"
    assert status["kept_run"] == "initial"
    assert status["rerun_learning_rate"] == "NA"
    assert (results / "cellbender_test.h5").read_text() == "fake raw for run 1\n"
    assert (results / "cellbender_test_report_initial.html").exists()
    assert not (results / "cellbender_test_report_rerun.html").exists()


def test_rerun_resolves_keeps_rerun(tmp_path, monkeypatch):
    _install_fake(tmp_path, monkeypatch, initial_phrase="rerun", rerun_phrase="normal")
    status, results = _run(tmp_path)

    assert status["outcome"] == "RERUN_RESOLVED"
    assert status["reran"] == "true"
    assert status["kept_run"] == "rerun"
    assert status["initial_learning_rate"] == "0.0001"
    assert status["rerun_learning_rate"] == "5e-05"  # exactly half
    # The kept output is the second (re-run) result, and both reports are archived.
    assert (results / "cellbender_test.h5").read_text() == "fake raw for run 2\n"
    assert (results / "cellbender_test_report_initial.html").exists()
    assert (results / "cellbender_test_report_rerun.html").exists()


def test_rerun_does_not_resolve_keeps_initial(tmp_path, monkeypatch):
    _install_fake(tmp_path, monkeypatch, initial_phrase="rerun", rerun_phrase="rerun")
    status, results = _run(tmp_path)

    assert status["outcome"] == "RERUN_DID_NOT_RESOLVE"
    assert status["reran"] == "true"
    assert status["kept_run"] == "initial"
    assert status["rerun_verdict"] == "rerun"
    # Kept output falls back to the initial run even though a re-run happened.
    assert (results / "cellbender_test.h5").read_text() == "fake raw for run 1\n"
    assert (results / "cellbender_test_report_rerun.html").exists()


def test_adaptive_off_does_not_rerun(tmp_path, monkeypatch):
    _install_fake(tmp_path, monkeypatch, initial_phrase="rerun", rerun_phrase="normal")
    status, results = _run(tmp_path, adaptive="false")

    assert status["outcome"] == "RERUN_SUGGESTED_BUT_ADAPTIVE_OFF"
    assert status["reran"] == "false"
    assert status["kept_run"] == "initial"
    assert (results / "cellbender_test.h5").read_text() == "fake raw for run 1\n"
    assert not (results / "cellbender_test_report_rerun.html").exists()


def test_unknown_initial_report_hard_errors(tmp_path, monkeypatch):
    _install_fake(tmp_path, monkeypatch, initial_phrase="unknown")
    with pytest.raises(SystemExit) as excinfo:
        _run(tmp_path)
    assert "could not classify" in str(excinfo.value)


def test_unknown_rerun_report_hard_errors(tmp_path, monkeypatch):
    _install_fake(tmp_path, monkeypatch, initial_phrase="rerun", rerun_phrase="unknown")
    with pytest.raises(SystemExit) as excinfo:
        _run(tmp_path)
    assert "re-run report" in str(excinfo.value)


# --- run-level summary ---------------------------------------------------------------------------

def _write_status(path, **fields):
    path.write_text("".join(f"{k}\t{v}\n" for k, v in fields.items()))
    return path


def test_summary_aggregates_status_files(tmp_path, capsys):
    a = _write_status(
        tmp_path / "a.txt", sample="alpha", outcome="NORMAL_FIRST_TRY",
        initial_learning_rate="0.0001", initial_verdict="normal", reran="false",
        rerun_learning_rate="NA", rerun_verdict="NA", kept_run="initial",
    )
    b = _write_status(
        tmp_path / "b.txt", sample="beta", outcome="RERUN_DID_NOT_RESOLVE",
        initial_learning_rate="0.0001", initial_verdict="rerun", reran="true",
        rerun_learning_rate="5e-05", rerun_verdict="rerun", kept_run="initial",
    )
    out = tmp_path / "summary.tsv"
    # Pass b before a to confirm the summary sorts rows by sample id.
    summary_mod.main(["--output", str(out), str(b), str(a)])

    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == list(summary_mod.COLUMNS)
    assert lines[1].startswith("alpha\t")
    assert lines[2].startswith("beta\t")

    printed = capsys.readouterr().out
    assert "2 sample(s)" in printed
    assert "RERUN_DID_NOT_RESOLVE: 1" in printed
    assert "beta" in printed  # named as re-ran and as unresolved
