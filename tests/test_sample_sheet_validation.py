import shutil
import subprocess
from pathlib import Path


REQUIRED_TENX_CHILDREN = [
    "filtered_feature_bc_matrix",
    "raw_feature_bc_matrix",
    "raw_feature_bc_matrix.h5",
]


def run_snakemake_with_sample_sheet(repo_root, sample_sheet, results_dir, *extra_args):
    snakemake = shutil.which("snakemake")
    assert snakemake is not None, "snakemake is not available on PATH"

    cmd = [
        snakemake,
        *extra_args,
        "--profile",
        "none",
        "--workflow-profile",
        "none",
        "--snakefile",
        "workflow/Snakefile",
        "--configfile",
        "config/config.yaml",
        "--config",
        f"sampleTable={sample_sheet}",
        f"resultsDir={results_dir}",
    ]
    return subprocess.run(
        cmd,
        cwd=repo_root,
        text=True,
        capture_output=True,
        check=False,
        timeout=120,
    )


def combined_output(result):
    return result.stdout + result.stderr


def write_sample_sheet(path, rows, header="sampleid\ttenx_datadir"):
    lines = [header]
    lines.extend(f"{sampleid}\t{tenx_datadir}" for sampleid, tenx_datadir in rows)
    path.write_text("\n".join(lines) + "\n")
    return path


def test_malformed_sample_sheet_missing_required_column_fails(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    sample_sheet = tmp_path / "malformed_samplesheet.tsv"
    sample_sheet.write_text("sampleid\nmissing_tenx_column\n")

    result = run_snakemake_with_sample_sheet(
        repo_root,
        sample_sheet,
        tmp_path / "results",
        "-np",
    )
    output = combined_output(result)

    assert result.returncode != 0
    assert "Invalid sample sheet" in output
    assert "missing required column(s): tenx_datadir" in output


def test_duplicate_sample_ids_fail_validation(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    sample_sheet = write_sample_sheet(
        tmp_path / "duplicate_samplesheet.tsv",
        [("duplicate", "testdata"), ("duplicate", "testdata")],
    )

    result = run_snakemake_with_sample_sheet(
        repo_root,
        sample_sheet,
        tmp_path / "results",
        "-np",
    )
    output = combined_output(result)

    assert result.returncode != 0
    assert "Invalid sample sheet" in output
    assert "duplicate sampleid value(s): duplicate" in output


def test_missing_tenx_datadir_fails_validation(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    missing_data_dir = tmp_path / "does_not_exist"
    sample_sheet = write_sample_sheet(
        tmp_path / "missing_data_samplesheet.tsv",
        [("missing_path", missing_data_dir)],
    )

    result = run_snakemake_with_sample_sheet(
        repo_root,
        sample_sheet,
        tmp_path / "results",
        "-np",
    )
    output = combined_output(result)

    assert result.returncode != 0
    assert "Invalid sample sheet" in output
    assert "tenx_datadir does not exist" in output
    assert str(missing_data_dir) in output


def test_incomplete_tenx_datadir_fails_with_missing_required_children(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    incomplete_data_dir = tmp_path / "incomplete_testdata"
    incomplete_data_dir.mkdir()
    (incomplete_data_dir / "filtered_feature_bc_matrix").mkdir()

    sample_sheet = write_sample_sheet(
        tmp_path / "incomplete_data_samplesheet.tsv",
        [("incomplete", incomplete_data_dir)],
    )

    result = run_snakemake_with_sample_sheet(
        repo_root,
        sample_sheet,
        tmp_path / "results",
        "-np",
    )
    output = combined_output(result)

    assert result.returncode != 0
    assert "Invalid sample sheet" in output
    assert "missing raw_feature_bc_matrix" in output
    assert "missing raw_feature_bc_matrix.h5" in output


def test_absolute_tenx_datadir_is_normalized_and_used_for_rule_inputs(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    absolute_testdata = (repo_root / "testdata").resolve()
    sample_sheet = write_sample_sheet(
        tmp_path / "absolute_data_samplesheet.tsv",
        [("absolute_path", absolute_testdata)],
    )

    result = run_snakemake_with_sample_sheet(
        repo_root,
        sample_sheet,
        tmp_path / "results",
        "-np",
    )
    output = combined_output(result)

    assert result.returncode == 0, output
    for child in REQUIRED_TENX_CHILDREN:
        assert str(absolute_testdata / child) in output


def test_sample_sheet_with_multiple_samples_expands_outputs_for_each_sample(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    sample_sheet = write_sample_sheet(
        tmp_path / "multi_sample_samplesheet.tsv",
        [("test_a", "testdata"), ("test_b", "testdata")],
    )
    results_dir = tmp_path / "results"

    result = run_snakemake_with_sample_sheet(
        repo_root,
        sample_sheet,
        results_dir,
        "-np",
    )
    output = combined_output(result)

    assert result.returncode == 0, output
    for sample_id in ["test_a", "test_b"]:
        assert f"wildcards: sample={sample_id}" in output
        assert str(results_dir / "seurat_filtered" / f"filtered_seurat_tenx_{sample_id}_markergenes.csv") in output
        assert str(results_dir / "cellbender" / f"cellbender_{sample_id}.h5") in output
        assert str(results_dir / "cellbender" / f"cellbender_{sample_id}_filtered.h5") in output


# --- mt_genes: an optional column naming a file of the sample's mitochondrial gene IDs ---------------------

MT_GENES_TEST = Path("testdata/mt_genes_test.txt")
# Scripts that build an object from the sample's features, and so take the sample's --mt-genes.
MT_GENES_SCRIPTS = ["tenx2seuratrds.R", "emptydrops.R", "cellbender2seurat.R", "soupx.R", "decontx.R"]
# Scripts that start from an upstream object and read the mapped features it records.
DOWNSTREAM_SCRIPTS = ["doubletfinder.R", "scdblfinder.R", "posthocfilter_threshold.R", "downsample_clusters.R"]


def write_mt_sample_sheet(path, rows):
    lines = ["sampleid\ttenx_datadir\tmt_genes"]
    lines.extend("\t".join(str(value) for value in row) for row in rows)
    path.write_text("\n".join(lines) + "\n")
    return path


def run_dry_run(repo_root, sample_sheet, results_dir, *config_entries):
    snakemake = shutil.which("snakemake")
    assert snakemake is not None, "snakemake is not available on PATH"
    cmd = [
        snakemake,
        "-np",
        "--profile", "none",
        "--workflow-profile", "none",
        "--snakefile", "workflow/Snakefile",
        "--configfile", "config/config.yaml",
        "--config",
        f"sampleTable={sample_sheet}",
        f"resultsDir={results_dir}",
        *config_entries,
    ]
    return subprocess.run(cmd, cwd=repo_root, text=True, capture_output=True, check=False, timeout=120)


def script_commands(output, script, sample=None):
    """The dry-run shell command lines (-p) that run workflow/scripts/<script> (for one sample, if given)."""
    commands = [line for line in output.splitlines() if f"Rscript workflow/scripts/{script} " in line]
    if sample is not None:
        commands = [line for line in commands if f"_{sample}" in line]
    return commands


def test_mt_genes_file_must_exist(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    missing = tmp_path / "no_such_list.txt"
    sample_sheet = write_mt_sample_sheet(tmp_path / "samplesheet.tsv", [("test", "testdata", missing)])

    result = run_dry_run(repo_root, sample_sheet, tmp_path / "results")
    output = combined_output(result)

    assert result.returncode != 0
    assert "Invalid sample sheet" in output
    assert f"sample test row 2: mt_genes file does not exist: {missing}" in output


def test_mt_genes_file_must_list_a_gene_id(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    empty_list = tmp_path / "comments_only.txt"
    empty_list.write_text("# no IDs here\n\n   \n#ENSMUSG00000025428\n")
    sample_sheet = write_mt_sample_sheet(tmp_path / "samplesheet.tsv", [("test", "testdata", empty_list)])

    result = run_dry_run(repo_root, sample_sheet, tmp_path / "results")
    output = combined_output(result)

    assert result.returncode != 0
    assert "Invalid sample sheet" in output
    assert "mt_genes file lists no gene IDs" in output


def test_mt_genes_list_must_name_one_of_the_sample_features(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    unknown = tmp_path / "unknown_ids.txt"
    unknown.write_text("ENSMUSG00000064341\nMT-ND1\n")
    sample_sheet = write_mt_sample_sheet(tmp_path / "samplesheet.tsv", [("test", "testdata", unknown)])

    result = run_dry_run(repo_root, sample_sheet, tmp_path / "results")
    output = combined_output(result)

    assert result.returncode != 0
    assert "[mt_genes] mitochondrial-gene problem(s)" in output
    assert "test: none of the 2 gene IDs in its mt_genes file" in output


def test_relative_mt_genes_path_reaches_the_scripts_that_build_objects(tmp_path):
    # Relative paths resolve against the working directory, as tenx_datadir's do.
    repo_root = Path(__file__).resolve().parents[1]
    sample_sheet = write_mt_sample_sheet(
        tmp_path / "samplesheet.tsv", [("test", "testdata", MT_GENES_TEST.as_posix())]
    )

    result = run_dry_run(
        repo_root, sample_sheet, tmp_path / "results", "posthoc_methods=[threshold]", "vars_to_regress=[percent.mt]"
    )
    output = combined_output(result)

    assert result.returncode == 0, output
    resolved = str((repo_root / MT_GENES_TEST).resolve())
    # Two of the list's seven IDs are not features of the test data: reported and skipped.
    assert "[mt_genes] WARNING: test: 2 of 7 gene IDs in its mt_genes file" in output
    assert "ENSMUSG00000064341, ENSMUSG00000064351" in output
    for script in MT_GENES_SCRIPTS:
        commands = script_commands(output, script)
        assert commands, f"no dry-run command runs {script}"
        for command in commands:
            assert f"--mt-genes={resolved} " in command, command
            assert "--vars-to-regress=percent.mt " in command, command
    for script in DOWNSTREAM_SCRIPTS:
        commands = script_commands(output, script)
        assert commands, f"no dry-run command runs {script}"
        for command in commands:
            assert "--mt-genes" not in command, command
            assert "--vars-to-regress=percent.mt " in command, command


def test_blank_mt_genes_keeps_the_pattern_for_that_sample(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    sample_sheet = write_mt_sample_sheet(
        tmp_path / "samplesheet.tsv",
        [("listed", "testdata", MT_GENES_TEST.as_posix()), ("unlisted", "testdata", "")],
    )

    result = run_dry_run(repo_root, sample_sheet, tmp_path / "results")
    output = combined_output(result)

    assert result.returncode == 0, output
    resolved = str((repo_root / MT_GENES_TEST).resolve())
    listed = script_commands(output, "tenx2seuratrds.R", sample="listed")
    unlisted = script_commands(output, "tenx2seuratrds.R", sample="unlisted")
    assert listed and all(f"--mt-genes={resolved} " in command for command in listed), listed
    assert unlisted and all("--mt-genes= " in command for command in unlisted), unlisted


def test_percent_mt_regression_ignores_an_excluded_sample_without_mitochondrial_features(tmp_path):
    repo_root = Path(__file__).resolve().parents[1]
    sample_sheet = write_mt_sample_sheet(
        tmp_path / "samplesheet.tsv",
        [("listed", "testdata", MT_GENES_TEST.as_posix()), ("unlisted", "testdata", "")],
    )

    excluded = run_dry_run(
        repo_root, sample_sheet, tmp_path / "results", "vars_to_regress=[percent.mt]", "excluded_samples=[unlisted]"
    )
    included = run_dry_run(repo_root, sample_sheet, tmp_path / "results", "vars_to_regress=[percent.mt]")

    assert excluded.returncode == 0, combined_output(excluded)
    assert included.returncode != 0
    assert "unlisted: percent.mt is in vars_to_regress, but no feature name in" in combined_output(included)
