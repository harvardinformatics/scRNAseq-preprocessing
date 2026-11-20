import os

import pytest


@pytest.fixture
def github_error():
    """Return a function that prints a failure as an error annotation when the tests run on GitHub Actions.

    The run's summary page shows annotations, and the API serves them without signing in, so a failure can
    be read there without opening the job log. pytest prints a failing test's output, annotation included.
    """

    def annotate(title, message, limit=6000):
        if os.environ.get("GITHUB_ACTIONS") != "true":
            return
        if len(message) > limit:
            message = message[:limit] + "\n... (cut short; the job log has the rest)"
        title = title.replace("%", "%25").replace(",", "%2C").replace(":", "%3A")
        message = message.replace("%", "%25").replace("\r", "%0D").replace("\n", "%0A")
        print(f"::error title={title}::{message}")

    return annotate


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
        "--run-downsample-scaling",
        action="store_true",
        default=False,
        help="run downsample_clusters on small and medium fixtures and check runtime scales reasonably",
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

