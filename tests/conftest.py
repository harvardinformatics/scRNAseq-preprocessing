def pytest_addoption(parser):
    parser.addoption(
        "--run-workflow",
        action="store_true",
        default=False,
        help="run the full Snakemake workflow on the test data",
    )
    parser.addoption(
        "--snakemake-conda-prefix",
        default=None,
        help="optional value to pass to Snakemake --conda-prefix for --run-workflow",
    )
    parser.addoption(
        "--run-downsample-rule",
        action="store_true",
        default=False,
        help="execute the downsample_clusters Snakemake rule on one downsampling test fixture",
    )
    parser.addoption(
        "--run-downsample-workflow",
        action="store_true",
        default=False,
        help="execute the full downsampling Snakemake workflow on downsampling testdata",
    )
    parser.addoption(
        "--run-conda-validation",
        action="store_true",
        default=False,
        help="create workflow conda environments and validate key package imports",
    )
    parser.addoption(
        "--conda-env-name",
        default=None,
        help="optional workflow/envs/*.yml filename to validate with --run-conda-validation",
    )
    parser.addoption(
        "--run-container-validation",
        action="store_true",
        default=False,
        help="pull workflow containers with docker, apptainer, or singularity",
    )
    parser.addoption(
        "--run-doubletfinder-install",
        action="store_true",
        default=False,
        help="execute the DoubletFinder GitHub install rule on test outputs",
    )

