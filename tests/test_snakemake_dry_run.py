import csv
import shutil
import subprocess
from pathlib import Path



SAMPLE_ID = "test"
TEST_SAMPLE_SHEET = Path("testdata/samplesheet_test.tsv")
RESULTS_DIR = "testdata/results"
DOWNSAMPLE_RESULTS_DIR = f"{RESULTS_DIR}/downsampling"
EXPECTED_OUTPUTS = Path(__file__).with_name("test_sample_rule_output_files.txt")


def run_snakemake(repo_root, *extra_args):
    sample_sheet = repo_root / TEST_SAMPLE_SHEET
    assert sample_sheet.exists(), f"missing test sample sheet: {sample_sheet}"

    snakemake = shutil.which("snakemake")
    assert snakemake is not None, "snakemake is not available on PATH"

    cmd = [
        snakemake,
        *extra_args,
        "--snakefile",
        "workflow/Snakefile",
        "--configfile",
        "config/config.yaml",
        "--config",
        f"sampleTable={TEST_SAMPLE_SHEET.as_posix()}",
        f"resultsDir={RESULTS_DIR}",
        f"downsampleResultsDir={DOWNSAMPLE_RESULTS_DIR}",
    ]
    return subprocess.run(
        cmd,
        cwd=repo_root,
        text=True,
        capture_output=True,
        check=False,
    )


def combined_output(result):
    return result.stdout + result.stderr


def test_testdata_dag_builds_with_dry_run():
    repo_root = Path(__file__).resolve().parents[1]
    result = run_snakemake(repo_root, "-np")
    output = combined_output(result)

    assert result.returncode == 0, output
    assert "Building DAG of jobs" in output
    assert "This was a dry-run" in output
    for rule_name in [
        "tenx2seuratrds",
        "emptydrops",
        "soupx",
        "soupx_emptydrops",
        "cellbender",
        "cellbender2seurat",
        "doubletfinder",
        "scdblfinder",
        "posthocfilter_threshold",
        "posthocfilter_mad",
        "combine_markers",
    ]:
        assert rule_name in output


def is_internal_workflow_output(output):
    return (
        "_markergenes_cluster" in output
        or output.endswith("_cluster_ids.txt")
        or output.endswith("_nclusters.txt")
        or output.endswith("_marker_manifest")
    )


def declared_sample_output_files(summary_text, expected_manifest):
    lines = summary_text.splitlines()
    header_index = next(
        i for i, line in enumerate(lines) if line.startswith("output_file\t")
    )
    reader = csv.DictReader(lines[header_index:], delimiter="\t")

    outputs = []
    for row in reader:
        output = row["output_file"]
        if f"_{SAMPLE_ID}" not in output and f"/{SAMPLE_ID}_" not in output:
            continue
        if row["status"] == "removed temp file":
            continue
        if is_internal_workflow_output(output):
            continue

        directory_children = [
            expected
            for expected in expected_manifest
            if expected.startswith(f"{output}/")
        ]
        if directory_children:
            outputs.extend(directory_children)
        else:
            outputs.append(output)

    return sorted(set(outputs))


def test_declared_test_sample_outputs_match_manifest():
    repo_root = Path(__file__).resolve().parents[1]
    expected = EXPECTED_OUTPUTS.read_text().splitlines()

    result = run_snakemake(repo_root, "--summary")
    output = combined_output(result)

    assert result.returncode == 0, output

    observed = declared_sample_output_files(output, expected)

    missing = sorted(set(expected) - set(observed))
    unexpected = sorted(set(observed) - set(expected))

    assert observed == expected, (
        "declared outputs differ from manifest\n"
        f"missing: {missing}\n"
        f"unexpected: {unexpected}"
    )
