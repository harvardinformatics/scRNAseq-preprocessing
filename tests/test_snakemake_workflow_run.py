import csv
import gzip
import hashlib
import math
import shutil
import subprocess
from pathlib import Path

import pytest


EXPECTED_OUTPUTS = Path(__file__).with_name("test_sample_rule_output_files.txt")
REFERENCE_OUTPUTS = Path(__file__).with_name("test_reference_output_files.txt")
REFERENCE_ROOT = Path(__file__).with_name("reference_outputs")
SEURAT_METADATA_COMPARATOR = Path(__file__).with_name("compare_seurat_metadata.R")
MARKER_ABS_TOLERANCE = 1e-8
MARKER_REL_TOLERANCE = 1e-6


def read_manifest(path):
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_seurat_object(path):
    return path.endswith(".rds")


def is_marker_table(path):
    return path.endswith("_markergenes.csv")


def is_emptydrops_matrix_component(path):
    return "/emptydrops/test_emptydrops_filtered_matrix/" in path and path.endswith(".gz")


def is_cellbender_h5(path):
    return path.startswith("testdata/results/cellbender/") and path.endswith(".h5")


def find_rscript_with_seurat(repo_root):
    candidates = []
    path_rscript = shutil.which("Rscript")
    if path_rscript:
        candidates.append(Path(path_rscript))
    candidates.extend(sorted((repo_root / ".snakemake" / "conda").glob("*/bin/Rscript")))

    checked = []
    for candidate in candidates:
        if candidate in checked or not candidate.exists():
            continue
        checked.append(candidate)
        result = subprocess.run(
            [
                str(candidate),
                "-e",
                "quit(status = ifelse(requireNamespace('Seurat', quietly = TRUE), 0, 1))",
            ],
            text=True,
            capture_output=True,
            check=False,
        )
        if result.returncode == 0:
            return candidate

    raise AssertionError("could not find an Rscript with the Seurat package installed")


def seurat_metadata_differences(repo_root, paths):
    """The comparator's report on every Seurat object whose metadata differs from its reference."""
    if not paths:
        return []

    rscript = find_rscript_with_seurat(repo_root)
    cmd = [str(rscript), str(SEURAT_METADATA_COMPARATOR)]
    for rel_path in paths:
        cmd.extend([rel_path, str(REFERENCE_ROOT / rel_path)])

    result = subprocess.run(cmd, cwd=repo_root, text=True, capture_output=True, check=False)
    if result.returncode == 0:
        return []
    return [(result.stdout + result.stderr).strip()]


def read_csv_table(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        return reader.fieldnames or [], rows


def numeric_value(value):
    if value == "":
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    if math.isnan(parsed):
        return None
    return parsed


def marker_value_matches(current, reference):
    current_numeric = numeric_value(current)
    reference_numeric = numeric_value(reference)
    if current_numeric is not None and reference_numeric is not None:
        return math.isclose(
            current_numeric,
            reference_numeric,
            rel_tol=MARKER_REL_TOLERANCE,
            abs_tol=MARKER_ABS_TOLERANCE,
        )
    return current == reference


def marker_rows_by_cluster_gene(rows):
    indexed = {}
    duplicates = set()
    for row in rows:
        key = (row.get("cluster"), row.get("genesymbol"))
        if key in indexed:
            duplicates.add(key)
        indexed[key] = row
    return indexed, duplicates


def marker_table_differences(repo_root, rel_path):
    """Every way one marker table differs from its reference: columns, (cluster, gene) rows and values."""
    current_columns, current_rows = read_csv_table(repo_root / rel_path)
    reference_columns, reference_rows = read_csv_table(REFERENCE_ROOT / rel_path)

    if current_columns != reference_columns:
        return [
            f"{rel_path}: marker columns differ\n"
            f"current: {current_columns}\nreference: {reference_columns}"
        ]
    missing_columns = [column for column in ("cluster", "genesymbol") if column not in current_columns]
    if missing_columns:
        return [f"{rel_path}: missing column(s) {', '.join(missing_columns)}"]

    current_index, current_duplicates = marker_rows_by_cluster_gene(current_rows)
    reference_index, reference_duplicates = marker_rows_by_cluster_gene(reference_rows)
    if current_duplicates or reference_duplicates:
        return [f"{rel_path}: duplicate cluster/genesymbol rows found"]

    differences = []
    missing = sorted(set(reference_index) - set(current_index))
    extra = sorted(set(current_index) - set(reference_index))
    if missing or extra:
        differences.append(
            f"{rel_path}: marker genes differ by cluster: {len(missing)} reference (cluster, gene) row(s) "
            f"missing, {len(extra)} extra; first missing {missing[:3]}, first extra {extra[:3]}"
        )

    shared = sorted(set(reference_index) & set(current_index))
    for column in current_columns:
        mismatched = [
            key
            for key in shared
            if not marker_value_matches(current_index[key][column], reference_index[key][column])
        ]
        if not mismatched:
            continue
        gaps = []
        for key in mismatched:
            current_numeric = numeric_value(current_index[key][column])
            reference_numeric = numeric_value(reference_index[key][column])
            if current_numeric is not None and reference_numeric is not None:
                gaps.append(abs(current_numeric - reference_numeric))
        largest = f"; largest abs diff {max(gaps):.6g}" if gaps else ""
        first = mismatched[0]
        differences.append(
            f"{rel_path}: marker column {column} differs in {len(mismatched)} of {len(shared)} shared rows"
            f"{largest}; first at {first}: "
            f"{current_index[first][column]!r} != {reference_index[first][column]!r}"
        )
    return differences


def read_gzip_contents(path):
    with gzip.open(path, "rb") as handle:
        return handle.read()


def emptydrops_matrix_differences(repo_root, paths):
    return [
        f"{rel_path}: decompressed matrix component differs"
        for rel_path in paths
        if read_gzip_contents(repo_root / rel_path) != read_gzip_contents(REFERENCE_ROOT / rel_path)
    ]


def h5_differences(repo_root, paths):
    differences = []
    for rel_path in paths:
        result = subprocess.run(
            ["h5diff", "-q", str(repo_root / rel_path), str(REFERENCE_ROOT / rel_path)],
            text=True,
            capture_output=True,
            check=False,
        )
        if result.returncode != 0:
            differences.append(
                (f"{rel_path}: HDF5 contents differ\n" + result.stdout + result.stderr).strip()
            )
    return differences


def byte_differences(repo_root, paths):
    return [
        f"{rel_path}: differs byte-for-byte from the reference"
        for rel_path in paths
        if file_sha256(repo_root / rel_path) != file_sha256(REFERENCE_ROOT / rel_path)
    ]


def reference_differences(repo_root):
    """Every way the current outputs differ from the reference snapshot (empty when they match).

    Every comparison runs to the end, so that one run of the full test-data workflow (which takes hours)
    shows all the differences instead of stopping at the first.
    """
    expected = read_manifest(REFERENCE_OUTPUTS)

    missing_current = [path for path in expected if not (repo_root / path).exists()]
    missing_reference = [path for path in expected if not (REFERENCE_ROOT / path).exists()]
    differences = [f"{path}: missing from the current outputs" for path in missing_current]
    differences += [f"{path}: missing from the reference outputs" for path in missing_reference]
    missing = set(missing_current) | set(missing_reference)
    present = [path for path in expected if path not in missing]

    seurat_objects = [path for path in present if is_seurat_object(path)]
    marker_tables = [path for path in present if is_marker_table(path)]
    matrix_components = [path for path in present if is_emptydrops_matrix_component(path)]
    h5_outputs = [path for path in present if is_cellbender_h5(path)]
    compared_otherwise = set(seurat_objects + marker_tables + matrix_components + h5_outputs)
    byte_identical = [path for path in present if path not in compared_otherwise]

    differences += seurat_metadata_differences(repo_root, seurat_objects)
    for path in marker_tables:
        differences += marker_table_differences(repo_root, path)
    differences += emptydrops_matrix_differences(repo_root, matrix_components)
    differences += h5_differences(repo_root, h5_outputs)
    differences += byte_differences(repo_root, byte_identical)
    return differences


def tail(text, lines=80):
    return "\n".join(text.splitlines()[-lines:])


def test_testdata_workflow_run(pytestconfig, github_error):
    if not pytestconfig.getoption("--run-workflow"):
        pytest.skip("use --run-workflow to run the full workflow on testdata")

    repo_root = Path(__file__).resolve().parents[1]
    cmd = ["bash", "tests/run_test_workflow.sh"]

    conda_prefix = pytestconfig.getoption("--snakemake-conda-prefix")
    if conda_prefix:
        cmd.extend(["--conda-prefix", conda_prefix])

    result = subprocess.run(
        cmd,
        cwd=repo_root,
        text=True,
        capture_output=True,
        check=False,
    )

    if result.returncode != 0:
        output = result.stdout + result.stderr
        github_error("Test-data workflow failed", tail(output))
        pytest.fail(output, pytrace=False)

    expected = read_manifest(EXPECTED_OUTPUTS)
    missing = [path for path in expected if not (repo_root / path).exists()]
    if missing:
        report = "missing expected workflow outputs:\n" + "\n".join(missing[:50])
        github_error("Test-data workflow outputs missing", report)
        pytest.fail(report, pytrace=False)

    differences = reference_differences(repo_root)
    if differences:
        report = "The test-data outputs differ from tests/reference_outputs:\n" + "\n".join(differences)
        github_error("Test-data outputs differ from the reference", report)
        pytest.fail(report, pytrace=False)
