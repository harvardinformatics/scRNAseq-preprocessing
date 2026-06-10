import csv
import os
import shutil
import subprocess
from pathlib import Path


TEST_SAMPLE_SHEET = Path("testdata/samplesheet_test.tsv")
CLUSTERS = ["0", "1", "3"]


FAKE_RSCRIPT = r"""#!/usr/bin/env bash
set -euo pipefail

script="${1:?missing script path}"
shift
log="${FAKE_RSCRIPT_LOG:?missing FAKE_RSCRIPT_LOG}"

case "$(basename "${script}")" in
    find_markers.R)
        rds="${1:?missing input rds}"
        cluster="${2:?missing cluster id}"
        output="${3:?missing output csv}"
        mkdir -p "$(dirname "${output}")"
        printf 'find_markers\t%s\t%s\t%s\n' "${cluster}" "${rds}" "${output}" >> "${log}"
        {
            printf 'genesymbol,p_val,avg_log2FC,pct.1,pct.2,p_val_adj,cluster,workflow\n'
            printf 'Gene%s,0.01,1.0,0.5,0.1,0.05,%s,checkpoint_test\n' "${cluster}" "${cluster}"
        } > "${output}"
        ;;
    combine_markers.R)
        output="${1:?missing output csv}"
        shift
        mkdir -p "$(dirname "${output}")"
        printf 'combine_markers\t%s\n' "$*" >> "${log}"
        : > "${output}"
        header_written=0
        for marker_csv in "$@"; do
            if [[ "${header_written}" -eq 0 ]]; then
                cat "${marker_csv}" >> "${output}"
                header_written=1
            else
                tail -n +2 "${marker_csv}" >> "${output}"
            fi
        done
        ;;
    *)
        echo "unexpected Rscript target: ${script}" >&2
        exit 2
        ;;
esac
"""


def write_fake_rscript(tmp_path):
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_rscript = fake_bin / "Rscript"
    fake_rscript.write_text(FAKE_RSCRIPT)
    fake_rscript.chmod(0o755)
    return fake_bin


def run_command(cmd, repo_root, env):
    return subprocess.run(
        cmd,
        cwd=repo_root,
        env=env,
        text=True,
        capture_output=True,
        check=False,
        timeout=120,
    )


def combined_output(result):
    return result.stdout + result.stderr


def read_marker_rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def test_marker_checkpoint_expands_find_markers_for_materialized_clusters(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    sample_sheet = repo_root / TEST_SAMPLE_SHEET
    assert sample_sheet.exists(), f"missing test sample sheet: {sample_sheet}"

    snakemake = shutil.which("snakemake")
    assert snakemake is not None, "snakemake is not available on PATH"

    fake_bin = write_fake_rscript(tmp_path)
    rscript_log = tmp_path / "fake_rscript.log"
    results_dir = tmp_path / "results"
    prefix = Path("checkpoint_marker") / "seurat_checkpoint_marker_test"
    prefix_path = results_dir / prefix
    marker_csv = prefix_path.with_name(prefix_path.name + "_markergenes.csv")
    cluster_ids = prefix_path.with_name(prefix_path.name + "_cluster_ids.txt")
    input_rds = prefix_path.with_suffix(".rds")
    manifest_dir = prefix_path.with_name(prefix_path.name + "_marker_manifest")

    prefix_path.parent.mkdir(parents=True)
    input_rds.write_text("dummy rds path consumed only by fake Rscript\n")
    cluster_ids.write_text("\n".join(CLUSTERS) + "\n")

    env = os.environ.copy()
    env["PATH"] = f"{fake_bin}{os.pathsep}{env.get('PATH', '')}"
    env["FAKE_RSCRIPT_LOG"] = str(rscript_log)

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
        "5",
        "--printshellcmds",
        "--snakefile",
        "workflow/Snakefile",
        "--configfile",
        "config/config.yaml",
        "--config",
        f"sampleTable={TEST_SAMPLE_SHEET.as_posix()}",
        f"resultsDir={results_dir.as_posix()}",
    ]

    result = run_command(cmd, repo_root, env)
    output = combined_output(result)
    assert result.returncode == 0, output
    assert "marker_manifest" in output
    assert "find_markers" in output
    assert "combine_markers" in output

    rows = read_marker_rows(marker_csv)
    observed_clusters = sorted(row["cluster"] for row in rows)
    assert observed_clusters == CLUSTERS
    assert [row["workflow"] for row in rows] == ["checkpoint_test"] * len(CLUSTERS)

    log_lines = rscript_log.read_text().splitlines()
    find_lines = [line for line in log_lines if line.startswith("find_markers\t")]
    combine_lines = [line for line in log_lines if line.startswith("combine_markers\t")]
    assert sorted(line.split("\t")[1] for line in find_lines) == CLUSTERS
    assert len(combine_lines) == 1
    for cluster in CLUSTERS:
        assert f"_markergenes_cluster{cluster}.csv" in combine_lines[0]

    assert not cluster_ids.exists(), "combine_markers should remove cluster id temp input"
    assert not manifest_dir.exists(), "combine_markers should remove checkpoint manifest directory"
