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


def assert_seurat_metadata_matches(repo_root, paths):
    if not paths:
        return

    rscript = find_rscript_with_seurat(repo_root)
    cmd = [str(rscript), str(SEURAT_METADATA_COMPARATOR)]
    for rel_path in paths:
        cmd.extend([str(repo_root / rel_path), str(REFERENCE_ROOT / rel_path)])

    result = subprocess.run(cmd, text=True, capture_output=True, check=False)
    assert result.returncode == 0, result.stdout + result.stderr


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


def assert_marker_value_matches(path, column, current, reference):
    current_numeric = numeric_value(current)
    reference_numeric = numeric_value(reference)
    if current_numeric is not None and reference_numeric is not None:
        assert math.isclose(
            current_numeric,
            reference_numeric,
            rel_tol=MARKER_REL_TOLERANCE,
            abs_tol=MARKER_ABS_TOLERANCE,
        ), (
            f"{path}: marker value differs for column {column}: "
            f"{current} != {reference}"
        )
        return

    assert current == reference, (
        f"{path}: marker value differs for column {column}: {current!r} != {reference!r}"
    )


def marker_rows_by_cluster_gene(path, rows):
    keys = []
    indexed = {}
    for row in rows:
        key = (row.get("cluster"), row.get("genesymbol"))
        keys.append(key)
        indexed[key] = row

    assert len(indexed) == len(rows), f"{path}: duplicate cluster/genesymbol rows found"
    return keys, indexed


def assert_marker_tables_match(repo_root, paths):
    for rel_path in paths:
        current_path = repo_root / rel_path
        reference_path = REFERENCE_ROOT / rel_path
        current_columns, current_rows = read_csv_table(current_path)
        reference_columns, reference_rows = read_csv_table(reference_path)

        assert current_columns == reference_columns, (
            f"{rel_path}: marker columns differ\n"
            f"current: {current_columns}\nreference: {reference_columns}"
        )
        assert "cluster" in current_columns, f"{rel_path}: missing cluster column"
        assert "genesymbol" in current_columns, f"{rel_path}: missing genesymbol column"

        current_keys, current_index = marker_rows_by_cluster_gene(rel_path, current_rows)
        reference_keys, reference_index = marker_rows_by_cluster_gene(rel_path, reference_rows)

        assert sorted(current_keys) == sorted(reference_keys), (
            f"{rel_path}: marker genes differ by cluster"
        )

        for key in sorted(reference_keys):
            current_row = current_index[key]
            reference_row = reference_index[key]
            for column in current_columns:
                assert_marker_value_matches(
                    rel_path,
                    column,
                    current_row[column],
                    reference_row[column],
                )


def read_gzip_contents(path):
    with gzip.open(path, "rb") as handle:
        return handle.read()


def assert_emptydrops_matrix_contents_match(repo_root, paths):
    for rel_path in paths:
        current = read_gzip_contents(repo_root / rel_path)
        reference = read_gzip_contents(REFERENCE_ROOT / rel_path)
        assert current == reference, f"{rel_path}: decompressed matrix component differs"


def assert_h5_contents_match(repo_root, paths):
    for rel_path in paths:
        result = subprocess.run(
            ["h5diff", "-q", str(repo_root / rel_path), str(REFERENCE_ROOT / rel_path)],
            text=True,
            capture_output=True,
            check=False,
        )
        assert result.returncode == 0, (
            f"{rel_path}: HDF5 contents differ\n" + result.stdout + result.stderr
        )


def assert_byte_identical(repo_root, paths):
    mismatches = []
    for rel_path in paths:
        current = repo_root / rel_path
        reference = REFERENCE_ROOT / rel_path
        if file_sha256(current) != file_sha256(reference):
            mismatches.append(rel_path)

    assert not mismatches, (
        "current outputs differ byte-for-byte from reference outputs:\n"
        + "\n".join(mismatches[:50])
    )


def assert_current_outputs_match_reference(repo_root):
    expected = read_manifest(REFERENCE_OUTPUTS)

    missing_current = [path for path in expected if not (repo_root / path).exists()]
    missing_reference = [
        path for path in expected if not (REFERENCE_ROOT / path).exists()
    ]
    assert not missing_current, (
        "missing current outputs listed in reference manifest:\n"
        + "\n".join(missing_current[:50])
    )
    assert not missing_reference, (
        "missing reference outputs listed in reference manifest:\n"
        + "\n".join(missing_reference[:50])
    )

    seurat_objects = [path for path in expected if is_seurat_object(path)]
    marker_tables = [path for path in expected if is_marker_table(path)]
    matrix_components = [
        path for path in expected if is_emptydrops_matrix_component(path)
    ]
    h5_outputs = [path for path in expected if is_cellbender_h5(path)]
    byte_identical = [
        path
        for path in expected
        if path not in set(seurat_objects + marker_tables + matrix_components + h5_outputs)
    ]

    assert_seurat_metadata_matches(repo_root, seurat_objects)
    assert_marker_tables_match(repo_root, marker_tables)
    assert_emptydrops_matrix_contents_match(repo_root, matrix_components)
    assert_h5_contents_match(repo_root, h5_outputs)
    assert_byte_identical(repo_root, byte_identical)


def test_testdata_workflow_run(pytestconfig):
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

    assert result.returncode == 0, result.stdout + result.stderr

    expected = read_manifest(EXPECTED_OUTPUTS)
    missing = [path for path in expected if not (repo_root / path).exists()]

    assert not missing, "missing expected workflow outputs:\n" + "\n".join(missing[:50])

    assert_current_outputs_match_reference(repo_root)
