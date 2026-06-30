# Tests

These tests are intended to run from the repository root in any environment where `snakemake` and `pytest` are available on `PATH`. How that environment is created or activated is site-specific; for example, an HPC may require loading a module before `conda activate` is available.

Example:

```bash
conda env create -f environment.yml
conda activate scrnaseq-preprocessing-tests
pytest tests
```

`python -m pytest tests` is also fine when `python` resolves to the same environment that provides `pytest` and `snakemake`. If it does not, the active shell is likely resolving commands from different environments. Confirm with:

```bash
which python
which pytest
which snakemake
```

The dry-run DAG test uses the small CellRanger-style input data in `testdata/` and the test-specific sample sheet `testdata/samplesheet_test.tsv`. The tests override `sampleTable` on the Snakemake command line, so they do not use the normal workflow `samplesheet.tsv`. They do not require `--conda-prefix`, a cluster profile, or any site-specific paths. Sample-sheet validation tests cover malformed headers, duplicate sample IDs, missing or incomplete 10x data paths, absolute `tenx_datadir` values, and multi-sample DAG expansion. Config validation tests cover required keys, allowed method values, threshold types/ranges, invalid-config failures, and dry-run DAG variants driven by method-list config changes. The default suite also runs `snakemake --lint` against the test-data configuration, which catches workflow-structure problems such as missing `log:` directives, helper functions embedded in rule files, long `run:` blocks that should live in scripts, path-composition warnings, missing rule-level conda/container declarations, and shell commands that directly interpolate global workflow variables instead of passing values through `params`.

GitHub Actions runs the default test suite via `.github/workflows/tests.yml`, using the same top-level `environment.yml` test runner environment and caching Snakemake-created rule environments under `.snakemake/conda`.

The default suite also includes the downsampling fixture/reference checks and a dry-run of `workflow_mode=downsample_only`. Optional downsampling execution tests can be run with `pytest tests/downsampling --run-downsample-rule -q` or `pytest tests/downsampling --run-downsample-workflow -q`, and both are available as manual GitHub Actions dispatch jobs.

The default test suite also includes a focused local rule-execution smoke test. It runs the real `tenx2seuratrds`, `find_markers`, and `combine_markers` rule chain against `testdata/`, using Snakemake's `--use-conda` support and writing outputs under pytest's temporary directory. This catches broken R package imports, script argument drift, invalid Seurat object creation, and marker CSV schema changes without submitting to SLURM.

The R output validator also checks that the Seurat object has at least 100 features and 100 cells; metadata rows match the cell count; barcode row names are present, unique, and nonempty; `orig.ident`, `nCount_RNA`, `nFeature_RNA`, `percent.mt`, and `seurat_clusters` metadata columns exist; RNA count and feature-count metadata values are finite and positive; mitochondrial percentages are finite and within `[0, 100]`; at least two clusters are present; PCA and UMAP reductions exist; the marker table is nonempty and has the expected columns; marker gene symbols are present and nonempty; marker numeric columns are finite; marker p-value and percent columns are within `[0, 1]`; marker clusters are present in the Seurat metadata; markers are reported for at least two clusters; and the marker `workflow` column matches the expected test workflow label. The test runner environment is defined in the repository-level `environment.yml`; the rule-specific R environment is still created by Snakemake from `workflow/envs/tenx2seuratrds.yml`. A separate lightweight checkpoint-expansion test uses a fake `Rscript` to materialize the `marker_manifest` checkpoint, verify that dynamic `find_markers` jobs are generated for each cluster id, and confirm that `combine_markers` receives the expected marker chunks.


## Optional conda and container validation

The default suite validates that workflow conda environment files are well formed, that rule-level conda references resolve to existing files, and that container declarations are recognizable. Expensive network-dependent validation is opt-in:

```bash
pytest tests --run-conda-validation
pytest tests --run-container-validation
pytest tests --run-doubletfinder-install
```

`--run-conda-validation` creates each `workflow/envs/*.yml` environment in a pytest temporary directory and checks key R/Python package imports. To validate only one environment, pass `--conda-env-name ENV_FILE`, for example `--conda-env-name soupx.yml`; GitHub Actions uses this selector to run conda validation as one matrix job per env. `--run-container-validation` pulls the CellBender container with Docker, Apptainer, or Singularity. `--run-doubletfinder-install` runs the networked Snakemake `install_doubletfinder` rule and confirms that `DoubletFinder` can be imported from the created rule environment.


## Optional full workflow run

The default tests build and inspect the DAG and run a focused local R-rule smoke test. To submit the full test-data workflow and verify all declared outputs against the reference snapshot, opt in explicitly:

```bash
pytest tests --run-workflow
```

The full-run test calls `tests/run_test_workflow.sh`, which uses `testdata/samplesheet_test.tsv` and overrides the workflow output directories with `resultsDir=testdata/results` and `downsampleResultsDir=testdata/results/downsampling`. The manifest in `tests/test_sample_rule_output_files.txt` is therefore written with paths under `testdata/results/`.

For testing, omit `--snakemake-conda-prefix` so Snakemake uses its default `.snakemake/conda` location under the repository root. The runner assumes that the current environment already provides `snakemake` on `PATH`.


## Reference outputs

The full workflow test compares regenerated files in `testdata/results/` against reference files under `tests/reference_outputs/`. The compared file list is in `tests/test_reference_output_files.txt`. Seurat `.rds` files are compared at the metadata-table level, marker CSVs are compared by columns and `(cluster, genesymbol)` rows with numeric tolerance, emptyDrops matrix files are compared after gzip decompression, CellBender H5 outputs are compared with `h5diff`, and remaining durable outputs are compared byte-for-byte.

To refresh the reference snapshot after intentionally changing workflow behavior, first run the full test workflow so `testdata/results/` contains the desired outputs, then run:

```bash
python tests/update_reference_outputs.py
```
