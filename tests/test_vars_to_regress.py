"""percent.mt and SCTransform's regression (config vars_to_regress, sample-sheet mt_genes), run in R.

The test data have no mitochondrial gene, so testdata/mt_genes_test.txt lists five nuclear-encoded
respiratory-chain genes in their place, plus two IDs that are not features (reported and skipped), and a copy
of the test matrix renames the same five genes mt-... or MT-..., for the "(?i)^mt-" pattern to find. One
Snakemake run builds the CellRanger object with the list and percent.mt regressed out (the rule passes the
options); the other CellRanger objects and the later steps run their scripts directly with Rscript, as the
workflow does. tests/check_percent_mt.R recomputes percent.mt independently and checks the regression. The
doubletfinder.R cases need an Rscript with DoubletFinder (first on PATH, or in .snakemake/conda once the
workflow's install_doubletfinder rule has installed it there) and are skipped without one.
"""

import gzip
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
MT_GENES_TEST = REPO_ROOT / "testdata" / "mt_genes_test.txt"
TEST_MATRIX = REPO_ROOT / "testdata" / "filtered_feature_bc_matrix"
CHECK_PERCENT_MT = Path(__file__).with_name("check_percent_mt.R")
CHECK_MAPPING = Path(__file__).with_name("check_mt_gene_mapping.R")
# Built with the "(?i)^mt-" pattern from the test data, which have no feature it matches: percent.mt is 0 in
# every cell.
NO_MT_OBJECT = (
    REPO_ROOT / "tests" / "reference_outputs" / "testdata" / "results" / "scdblfinder"
    / "seurat_scdblfinder_soupx_tenx_test.rds"
)
SCRIPTS = Path("workflow/scripts")
# Sets percent.mt as the doublet-removal scripts computed it before every step computed it on the RNA counts:
# on the object's SCT assay, its default.
SET_SCT_BASED_PERCENT_MT = """\
suppressPackageStartupMessages(library(Seurat))
args <- commandArgs(trailingOnly = TRUE)
x <- readRDS(args[1])
listed <- x@misc$percent_mt_features
x$percent.mt <- if (is.null(listed)) {
  PercentageFeatureSet(x, pattern = "(?i)^mt-", assay = "SCT")
} else {
  PercentageFeatureSet(x, features = intersect(listed, rownames(x[["SCT"]])), assay = "SCT")
}
saveRDS(x, args[2])
"""


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


def find_rscript_with(*packages):
    """The first Rscript on PATH, else in the R environments Snakemake created, that has every package; or None."""
    candidates = []
    path_rscript = shutil.which("Rscript")
    if path_rscript:
        candidates.append(Path(path_rscript))
    candidates.extend(sorted((REPO_ROOT / ".snakemake" / "conda").glob("*/bin/Rscript")))
    installed = " && ".join(f"requireNamespace('{package}', quietly = TRUE)" for package in packages)
    for candidate in dict.fromkeys(candidates):
        if not candidate.exists():
            continue
        probe = run([candidate, "-e", f"quit(status = ifelse({installed}, 0, 1))"])
        if probe.returncode == 0:
            return candidate
    return None


def find_rscript_with_seurat():
    rscript = find_rscript_with("Seurat")
    assert rscript is not None, "could not find an Rscript with the Seurat package installed"
    return rscript


def strip_r_comments(code):
    """R code with its comments (from a # outside a string or backquoted name to the end of the line) removed."""
    kept = []
    quote = None
    in_comment = False
    i = 0
    while i < len(code):
        char = code[i]
        if in_comment:
            if char == "\n":
                in_comment = False
                kept.append(char)
        elif quote:
            kept.append(char)
            if char == "\\":
                kept.append(code[i + 1:i + 2])
                i += 1
            elif char == quote:
                quote = None
        elif char in "\"'`":
            quote = char
            kept.append(char)
        elif char == "#":
            in_comment = True
        else:
            kept.append(char)
        i += 1
    return "".join(kept)


def closing_index(code, start):
    """The index just past the bracket that closes the one at code[start], in R code without comments."""
    opening = code[start]
    closing = {"(": ")", "{": "}", "[": "]"}[opening]
    depth = 0
    quote = None
    i = start
    while i < len(code):
        char = code[i]
        if quote:
            if char == "\\":
                i += 1
            elif char == quote:
                quote = None
        elif char in "\"'`":
            quote = char
        elif char == opening:
            depth += 1
        elif char == closing:
            depth -= 1
            if depth == 0:
                return i + 1
        i += 1
    raise AssertionError(f"unbalanced {opening!r} at offset {start}")


def listed_ids():
    lines = (line.strip() for line in MT_GENES_TEST.read_text().splitlines())
    return {line for line in lines if line and not line.startswith("#")}


def write_gz(path, text):
    with gzip.open(path, "wt") as handle:
        handle.write(text)


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


@pytest.fixture(scope="module")
def pattern_matrix(regressed_build):
    """The test matrix with the five listed genes renamed so that the "(?i)^mt-" pattern finds them."""
    matrix_dir = regressed_build["tmp"] / "pattern_matrix"
    matrix_dir.mkdir()
    for name in ("barcodes.tsv.gz", "matrix.mtx.gz"):
        shutil.copyfile(TEST_MATRIX / name, matrix_dir / name)
    ids = listed_ids()
    lines = []
    renamed = 0
    with gzip.open(TEST_MATRIX / "features.tsv.gz", "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] in ids:
                fields[1] = ("mt-" if renamed % 2 == 0 else "MT-") + fields[1]
                renamed += 1
            lines.append("\t".join(fields) + "\n")
    assert renamed == 5
    write_gz(matrix_dir / "features.tsv.gz", "".join(lines))
    return matrix_dir


@pytest.fixture(scope="module")
def sct_objects(regressed_build, pattern_matrix):
    """CellRanger objects built without regression, which default to their SCT assay as every saved object does:
    one whose percent.mt counts the genes the pattern finds, one whose percent.mt counts the same genes, listed."""
    tmp = regressed_build["tmp"]
    objects = {}
    for name, matrix, mt_genes in (("pattern", pattern_matrix, ""), ("list", TEST_MATRIX, MT_GENES_TEST)):
        rds = tmp / f"tenx_{name}.rds"
        result = run_script(
            regressed_build, "tenx2seuratrds.R", matrix, rds, tmp / f"tenx_{name}_n.txt", tmp / f"tenx_{name}_ids.txt",
            f"--mt-genes={mt_genes}", "--vars-to-regress=", "--sample=test",
        )
        assert result.returncode == 0, output_of(result)
        objects[name] = {"rds": rds, "matrix": matrix, "log": output_of(result)}
    return objects


@pytest.fixture(scope="module")
def doubletfinder_rscript():
    """An Rscript with DoubletFinder, or None. The R environment Snakemake creates has none until the workflow's
    install_doubletfinder rule installs it from GitHub."""
    return find_rscript_with("Seurat", "DoubletFinder")


@pytest.fixture(scope="module")
def doublet_step(regressed_build, sct_objects):
    """Runs a doublet-removal script once per (script, CellRanger object, --vars-to-regress)."""
    outputs = {}

    def run_step(script, build, regressed, rscript=None):
        key = (script, build, regressed)
        if key not in outputs:
            tmp = regressed_build["tmp"]
            stem = f"{Path(script).stem}_{build}_{regressed or 'none'}"
            rds = tmp / f"{stem}.rds"
            result = run([
                rscript or regressed_build["rscript"], SCRIPTS / script, sct_objects[build]["rds"], rds,
                tmp / f"{stem}_n.txt", tmp / f"{stem}_ids.txt", f"--vars-to-regress={regressed or ''}", "--sample=test",
            ])
            outputs[key] = {"rds": rds, "result": result}
        return outputs[key]

    return run_step


def test_only_add_percent_mt_computes_percent_mt_and_on_the_rna_assay():
    # Every percent.mt a script stores comes from add_percent_mt() in regression_utils.R, on the RNA assay. A call
    # of PercentageFeatureSet() without assay = "RNA" computes it on the object's default assay, which is SCT in
    # the doublet-removal steps' inputs.
    scripts = REPO_ROOT / "workflow" / "scripts"
    utils = scripts / "regression_utils.R"
    utils_code = strip_r_comments(utils.read_text())
    definition = re.search(r"^add_percent_mt\s*<-\s*function\s*\(", utils_code, re.MULTILINE)
    assert definition, "regression_utils.R does not define add_percent_mt()"
    body_start = utils_code.index("{", definition.end())
    body_end = closing_index(utils_code, body_start)

    def where(path, code, offset):
        return f"{path.relative_to(REPO_ROOT)}:{code.count(chr(10), 0, offset) + 1}"

    def in_add_percent_mt(path, offset):
        return path == utils and body_start < offset < body_end

    calls = []
    for path in sorted(scripts.rglob("*.R")):
        code = strip_r_comments(path.read_text())
        for call in re.finditer(r"\bPercentageFeatureSet\s*\(", code):
            site = where(path, code, call.start())
            assert in_add_percent_mt(path, call.start()), f"{site} calls PercentageFeatureSet() outside add_percent_mt()"
            arguments = code[call.end() - 1:closing_index(code, call.end() - 1)]
            assert re.search(r'\bassay\s*=\s*"RNA"', arguments), f'{site}: PercentageFeatureSet{arguments} lacks assay = "RNA"'
            calls.append(site)
        for store in re.finditer(r'\[\[\s*"percent\.mt"\s*\]\]\s*<-|\$\s*percent\.mt\s*<-|col\.name\s*=\s*"percent\.mt"', code):
            site = where(path, code, store.start())
            assert in_add_percent_mt(path, store.start()), f"{site} stores percent.mt outside add_percent_mt()"
    assert calls, "add_percent_mt() does not call PercentageFeatureSet()"


def test_rule_maps_the_list_and_regresses_percent_mt(regressed_build):
    log = regressed_build["log"]
    assert "tenx2seuratrds (sample test): 2 of 7 mt_genes IDs are not in column 1 of" in log
    assert "skipped: ENSMUSG00000064341, ENSMUSG00000064351" in log
    assert "tenx2seuratrds (sample test): 5 mt_genes ID(s) from" in log
    assert "tenx2seuratrds (sample test): SCTransform regresses out percent.mt" in log
    check(regressed_build, "built", regressed_build["rds"], TEST_MATRIX, MT_GENES_TEST, "percent.mt")


def test_without_regression_percent_mt_counts_the_pattern_or_the_list(regressed_build, sct_objects):
    for name, mt_genes in (("pattern", "pattern"), ("list", MT_GENES_TEST)):
        built = sct_objects[name]
        assert "SCTransform regresses out" not in built["log"]
        check(regressed_build, "built", built["rds"], built["matrix"], mt_genes, "none")


DOUBLET_CASES = [
    (script, build, regressed)
    for script in ("scdblfinder.R", "doubletfinder.R")
    for build in ("pattern", "list")
    for regressed in (None, "percent.mt")
]


@pytest.mark.parametrize(
    "script, build, regressed",
    DOUBLET_CASES,
    ids=[f"{Path(script).stem}-{build}{'-regressed' if regressed else ''}" for script, build, regressed in DOUBLET_CASES],
)
def test_doublet_removal_computes_percent_mt_on_the_rna_counts(
    regressed_build, sct_objects, doublet_step, doubletfinder_rscript, script, build, regressed
):
    # The input defaults to its SCT assay; percent.mt is computed on the RNA counts all the same, with or
    # without regression, from the features the input records or with the pattern.
    rscript = None
    if script == "doubletfinder.R":
        if doubletfinder_rscript is None:
            pytest.skip(
                "no Rscript with Seurat and DoubletFinder on PATH or in .snakemake/conda: the R environment Snakemake "
                "creates lacks DoubletFinder until the workflow's install_doubletfinder rule installs it; put an "
                "Rscript that has it first on PATH"
            )
        rscript = doubletfinder_rscript
    step = doublet_step(script, build, regressed, rscript)
    output = output_of(step["result"])
    assert step["result"].returncode == 0, output
    context = Path(script).stem
    if build == "list":
        assert f"{context} (sample test): percent.mt from 5 of the sample's 5 mt_genes feature(s) (RNA assay)" in output
    elif regressed:
        assert f'{context} (sample test): percent.mt from the 5 feature(s) matching "(?i)^mt-" (RNA assay)' in output
    assert (f"{context} (sample test): SCTransform regresses out percent.mt" in output) == bool(regressed)
    checked = check(regressed_build, "carried", step["rds"], sct_objects[build]["rds"], regressed or "none")
    assert "upstream default assay: SCT;" in checked


@pytest.mark.parametrize("build", ["pattern", "list"])
def test_threshold_filter_recomputes_percent_mt_on_the_rna_counts(regressed_build, doublet_step, build):
    # From an upstream object whose percent.mt was computed on its SCT assay, as the doublet-removal steps
    # computed it before: the threshold filter must not pass that value on.
    tmp = regressed_build["tmp"]
    singlets = doublet_step("scdblfinder.R", build, None)
    assert singlets["result"].returncode == 0, output_of(singlets["result"])
    set_sct_based = tmp / "set_sct_based_percent_mt.R"
    set_sct_based.write_text(SET_SCT_BASED_PERCENT_MT)
    sct_based = tmp / f"scdblfinder_{build}_sct_based.rds"
    result = run([regressed_build["rscript"], set_sct_based, singlets["rds"], sct_based])
    assert result.returncode == 0, output_of(result)

    filtered = tmp / f"posthocfilter_threshold_{build}.rds"
    result = run_script(
        regressed_build, "posthocfilter_threshold.R", sct_based, filtered, "1", "1", tmp / f"pt_{build}_n.txt",
        tmp / f"pt_{build}_ids.txt", "--vars-to-regress=", "--sample=test",
    )
    assert result.returncode == 0, output_of(result)
    check(regressed_build, "carried", filtered, sct_based, "none")


def test_threshold_filter_and_downsampling_regress_an_rna_percent_mt(regressed_build, doublet_step):
    tmp = regressed_build["tmp"]
    singlets = doublet_step("scdblfinder.R", "list", "percent.mt")
    assert singlets["result"].returncode == 0, output_of(singlets["result"])

    filtered = tmp / "posthocfilter_threshold_regressed.rds"
    result = run_script(
        regressed_build, "posthocfilter_threshold.R", singlets["rds"], filtered, "1", "1", tmp / "ptr_n.txt",
        tmp / "ptr_ids.txt", "--vars-to-regress=percent.mt", "--sample=test",
    )
    assert result.returncode == 0, output_of(result)
    assert "posthocfilter_threshold (sample test): percent.mt from 5 of the sample's 5 mt_genes feature(s) (RNA assay)" in output_of(result)
    check(regressed_build, "carried", filtered, singlets["rds"], "percent.mt")

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

    # So does a step whose input defaults to its SCT assay: the check, like percent.mt, is on the RNA assay.
    rds = tmp / "no_mt_scdblfinder.rds"
    result = run_script(
        regressed_build, "scdblfinder.R", NO_MT_OBJECT, rds, tmp / "no_mt_sc_n.txt", tmp / "no_mt_sc_ids.txt",
        "--vars-to-regress=percent.mt", "--sample=test",
    )
    assert result.returncode != 0
    assert (
        'scdblfinder (sample test): percent.mt is in vars_to_regress, but no feature of the RNA assay matches'
    ) in output_of(result)
    assert not rds.exists()

    # And the threshold filter, which recomputes it.
    rds = tmp / "no_mt_posthoc.rds"
    result = run_script(
        regressed_build, "posthocfilter_threshold.R", NO_MT_OBJECT, rds, "1", "1", tmp / "no_mt_posthoc_n.txt",
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
