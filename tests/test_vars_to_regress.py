"""percent.mt from an mt_genes list and SCTransform's regression (config vars_to_regress), run in R.

The test data have no mitochondrial gene, so testdata/mt_genes_test.txt lists five nuclear-encoded
respiratory-chain genes in their place, plus two IDs that are not features (reported and skipped). One
Snakemake run builds the CellRanger object with that list and percent.mt regressed out (the rule passes the
options); the later steps then run their scripts directly with Rscript from that object, as the workflow does.
tests/check_percent_mt.R recomputes percent.mt independently and checks the regression.
"""

import gzip
import os
import shutil
import subprocess
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
MT_GENES_TEST = REPO_ROOT / "testdata" / "mt_genes_test.txt"
TEST_MATRIX = REPO_ROOT / "testdata" / "filtered_feature_bc_matrix"
CHECK_PERCENT_MT = Path(__file__).with_name("check_percent_mt.R")
CHECK_MAPPING = Path(__file__).with_name("check_mt_gene_mapping.R")
# Built with the "(?i)^mt-" pattern, which finds nothing in the test data: percent.mt is 0 in every cell.
PATTERN_OBJECT = (
    REPO_ROOT / "tests" / "reference_outputs" / "testdata" / "results" / "scdblfinder"
    / "seurat_scdblfinder_soupx_tenx_test.rds"
)
SCRIPTS = Path("workflow/scripts")


def r_env():
    env = {key: value for key, value in os.environ.items() if key not in {"R_LIBS", "R_LIBS_USER", "R_LIBS_SITE"}}
    env["SCRNASEQ_PREPROCESS_SEED"] = "12345"
    env["SCRNASEQ_DOWNSAMPLE_SEED"] = "12345"
    return env


def run(cmd, timeout=1200):
    return subprocess.run(
        [str(part) for part in cmd], cwd=REPO_ROOT, env=r_env(), text=True, capture_output=True, check=False,
        timeout=timeout,
    )


def output_of(result):
    return result.stdout + result.stderr


def find_rscript_with_seurat():
    candidates = []
    path_rscript = shutil.which("Rscript")
    if path_rscript:
        candidates.append(Path(path_rscript))
    candidates.extend(sorted((REPO_ROOT / ".snakemake" / "conda").glob("*/bin/Rscript")))
    for candidate in dict.fromkeys(candidates):
        if not candidate.exists():
            continue
        probe = run([candidate, "-e", "quit(status = ifelse(requireNamespace('Seurat', quietly = TRUE), 0, 1))"])
        if probe.returncode == 0:
            return candidate
    raise AssertionError("could not find an Rscript with the Seurat package installed")


@pytest.fixture(scope="module")
def regressed_build(tmp_path_factory):
    """The tenx2seuratrds rule's object for the test sample, with the test list and percent.mt regressed out."""
    snakemake = shutil.which("snakemake")
    assert snakemake is not None, "snakemake is not available on PATH"
    tmp = tmp_path_factory.mktemp("vars_to_regress")
    sample_sheet = tmp / "samplesheet.tsv"
    sample_sheet.write_text("sampleid\ttenx_datadir\tmt_genes\ntest\ttestdata\ttestdata/mt_genes_test.txt\n")
    results = tmp / "results"
    rds = results / "seurat_filtered" / "filtered_seurat_tenx_test.rds"
    cmd = [
        snakemake, rds,
        "--profile", "none", "--workflow-profile", "none", "--executor", "local",
        "--cores", "1", "--jobs", "1", "--latency-wait", "30", "--use-conda", "--show-failed-logs",
        "--snakefile", "workflow/Snakefile", "--configfile", "config/config.yaml",
        "--config", f"sampleTable={sample_sheet}", f"resultsDir={results}", "vars_to_regress=[percent.mt]",
    ]
    conda_frontend = os.environ.get("SNAKEMAKE_CONDA_FRONTEND")
    if conda_frontend is None and shutil.which("mamba"):
        conda_frontend = "mamba"
    if conda_frontend:
        cmd.extend(["--conda-frontend", conda_frontend])
    result = run(cmd, timeout=3600)
    assert result.returncode == 0, output_of(result)
    assert rds.exists(), f"missing Seurat RDS: {rds}"
    log = (results / "logs" / "tenx2seuratrds" / "tenx2seuratrds_test.log").read_text()
    return {"rds": rds, "log": log, "tmp": tmp, "rscript": find_rscript_with_seurat()}


def run_script(build, script, *args):
    return run([build["rscript"], SCRIPTS / script, *args])


def check(build, *args):
    result = run([build["rscript"], CHECK_PERCENT_MT, *args])
    assert result.returncode == 0, output_of(result)
    return output_of(result)


def test_rule_maps_the_list_and_regresses_percent_mt(regressed_build):
    log = regressed_build["log"]
    assert "tenx2seuratrds (sample test): 2 of 7 mt_genes IDs are not in column 1 of" in log
    assert "skipped: ENSMUSG00000064341, ENSMUSG00000064351" in log
    assert "tenx2seuratrds (sample test): 5 mt_genes ID(s) from" in log
    assert "tenx2seuratrds (sample test): SCTransform regresses out percent.mt" in log
    check(regressed_build, "built", regressed_build["rds"], TEST_MATRIX, MT_GENES_TEST, "percent.mt")


def test_list_without_regression_only_changes_percent_mt(regressed_build):
    tmp = regressed_build["tmp"]
    rds = tmp / "list_only.rds"
    result = run_script(
        regressed_build, "tenx2seuratrds.R", TEST_MATRIX, rds, tmp / "list_only_n.txt", tmp / "list_only_ids.txt",
        f"--mt-genes={MT_GENES_TEST}", "--vars-to-regress=", "--sample=test",
    )
    assert result.returncode == 0, output_of(result)
    assert "SCTransform regresses out" not in output_of(result)
    check(regressed_build, "built", rds, TEST_MATRIX, MT_GENES_TEST, "none")


def test_later_steps_regress_an_rna_percent_mt_from_the_recorded_features(regressed_build):
    tmp = regressed_build["tmp"]
    upstream = regressed_build["rds"]

    # scDblFinder's input defaults to its SCT assay, as DoubletFinder's does; the percent.mt it regresses out
    # is computed on the RNA counts, from the features its input records.
    singlets = tmp / "scdblfinder.rds"
    result = run_script(
        regressed_build, "scdblfinder.R", upstream, singlets, tmp / "sc_n.txt", tmp / "sc_ids.txt",
        "--vars-to-regress=percent.mt", "--sample=test",
    )
    assert result.returncode == 0, output_of(result)
    assert "scdblfinder (sample test): percent.mt from 5 of the sample's 5 mt_genes feature(s) (RNA assay)" in output_of(result)
    checked = check(regressed_build, "carried", singlets, upstream, "rna", "percent.mt")
    assert "upstream default assay: SCT;" in checked

    # The threshold filter recomputes the percent.mt it regresses out on the RNA counts too.
    filtered = tmp / "posthocfilter_threshold.rds"
    result = run_script(
        regressed_build, "posthocfilter_threshold.R", singlets, filtered, "1", "1", tmp / "pt_n.txt",
        tmp / "pt_ids.txt", "--vars-to-regress=percent.mt", "--sample=test",
    )
    assert result.returncode == 0, output_of(result)
    assert "posthocfilter_threshold (sample test): percent.mt from 5 of the sample's 5 mt_genes feature(s) (RNA assay)" in output_of(result)
    check(regressed_build, "carried", filtered, singlets, "rna", "percent.mt")

    # Every downsampling replicate regresses it out too.
    table = tmp / "downsampling.tsv"
    result = run_script(
        regressed_build, "downsample_clusters.R", filtered, table, "1", "0.8", "--vars-to-regress=percent.mt",
        "--sample=test",
    )
    assert result.returncode == 0, output_of(result)
    assert "downsample_clusters (sample test): percent.mt from 5 of the sample's 5 mt_genes feature(s) (RNA assay)" in output_of(result)
    assert "downsample_clusters (sample test): SCTransform regresses out percent.mt" in output_of(result)
    assert table.read_text().startswith("clusterid\tmax_jaccard\tbootstrap_number\n")


def test_percent_mt_regression_without_a_mitochondrial_feature_fails(regressed_build):
    tmp = regressed_build["tmp"]
    rds = tmp / "no_mt.rds"
    result = run_script(
        regressed_build, "tenx2seuratrds.R", TEST_MATRIX, rds, tmp / "no_mt_n.txt", tmp / "no_mt_ids.txt",
        "--vars-to-regress=percent.mt", "--sample=test",
    )
    output = output_of(result)
    assert result.returncode != 0
    assert (
        'tenx2seuratrds (sample test): percent.mt is in vars_to_regress, but no feature of the RNA assay '
        'matches "(?i)^mt-"'
    ) in output
    assert "sample sheet's mt_genes column" in output
    assert not rds.exists()

    # So does a step whose input defaults to its SCT assay: the check, like the percent.mt it guards, is on RNA.
    rds = tmp / "no_mt_scdblfinder.rds"
    result = run_script(
        regressed_build, "scdblfinder.R", PATTERN_OBJECT, rds, tmp / "no_mt_sc_n.txt", tmp / "no_mt_sc_ids.txt",
        "--vars-to-regress=percent.mt", "--sample=test",
    )
    assert result.returncode != 0
    assert (
        'scdblfinder (sample test): percent.mt is in vars_to_regress, but no feature of the RNA assay matches'
    ) in output_of(result)
    assert not rds.exists()

    # And a step that keeps its input's percent.mt, computed with the pattern.
    rds = tmp / "no_mt_posthoc.rds"
    result = run_script(
        regressed_build, "posthocfilter_threshold.R", PATTERN_OBJECT, rds, "1", "1", tmp / "no_mt_posthoc_n.txt",
        tmp / "no_mt_posthoc_ids.txt", "--vars-to-regress=percent.mt", "--sample=test",
    )
    assert result.returncode != 0
    assert "posthocfilter_threshold (sample test): percent.mt is in vars_to_regress, but no feature" in output_of(result)
    assert not rds.exists()


def test_list_naming_no_feature_and_unknown_variables_fail(regressed_build):
    tmp = regressed_build["tmp"]
    unknown = tmp / "unknown_ids.txt"
    unknown.write_text("ENSMUSG00000064341\nMT-ND1\n")
    result = run_script(
        regressed_build, "tenx2seuratrds.R", TEST_MATRIX, tmp / "unknown.rds", tmp / "unknown_n.txt",
        tmp / "unknown_ids_out.txt", f"--mt-genes={unknown}", "--sample=test",
    )
    assert result.returncode != 0
    assert "tenx2seuratrds (sample test): none of the 2 gene IDs in the mt_genes file" in output_of(result)

    result = run_script(
        regressed_build, "tenx2seuratrds.R", TEST_MATRIX, tmp / "bad_var.rds", tmp / "bad_var_n.txt",
        tmp / "bad_var_ids.txt", "--vars-to-regress=percent.mito",
    )
    assert result.returncode != 0
    assert "--vars-to-regress: unknown variable(s) percent.mito" in output_of(result)


def write_gz(path, text):
    with gzip.open(path, "wt") as handle:
        handle.write(text)


def test_list_maps_to_the_row_names_seurat_gives(regressed_build, tmp_path):
    # Duplicated names get make.unique suffixes (Read10X), and "_" and "|" become "-" (CreateSeuratObject).
    matrix_dir = tmp_path / "matrix"
    matrix_dir.mkdir()
    features = [
        ("ID1", "Dup"), ("ID2", "Dup"), ("ID3", "mt_x"), ("ID4", "a|b"), ("ID5", "Dup"), ("ID6", "Plain"), ("ID7", "Other"),
    ]
    write_gz(matrix_dir / "features.tsv.gz", "".join(f"{i}\t{n}\tGene Expression\n" for i, n in features))
    write_gz(matrix_dir / "barcodes.tsv.gz", "AAAC-1\nAAAG-1\nAAAT-1\n")
    entries = [(gene, cell, gene + cell) for gene in range(1, 8) for cell in range(1, 4)]
    write_gz(
        matrix_dir / "matrix.mtx.gz",
        "%%MatrixMarket matrix coordinate integer general\n"
        f"7 3 {len(entries)}\n" + "".join(f"{g} {c} {v}\n" for g, c, v in entries),
    )
    mt_genes = tmp_path / "mt_genes.txt"
    mt_genes.write_text("# comment\nID2\n\n ID3 \nID4\r\nID5\nNOT_A_FEATURE\n")

    result = run([regressed_build["rscript"], CHECK_MAPPING, matrix_dir, mt_genes, "Dup.1", "mt-x", "a-b", "Dup.2"])
    assert result.returncode == 0, output_of(result)
    assert "1 of 5 mt_genes IDs are not in column 1" in output_of(result)
