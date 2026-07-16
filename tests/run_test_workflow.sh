#!/usr/bin/env bash
#SBATCH -p sapphire,shared
#SBATCH -e testdata_workflow_%A.err
#SBATCH -o testdata_workflow_%A.out
#SBATCH -J testdata_workflow
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem 8000
#SBATCH -t 72:00:00

set -euo pipefail

usage() {
    echo "Usage: $0 [--conda-prefix PATH] [additional snakemake args...]"
    echo "Set TEST_WORKFLOW_SNAKEMAKE_ARGS to append whitespace-delimited Snakemake args before command-line extras."
    echo "Set TEST_WORKFLOW_PROFILE_ARGS to replace the default Slurm/cannon profile args."
}

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "${script_dir}/.." && pwd)
cd "${repo_root}"

snakemake_conda_prefix=""
extra_snakemake_args=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --conda-prefix)
            if [[ $# -lt 2 ]]; then
                usage >&2
                exit 2
            fi
            snakemake_conda_prefix="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            extra_snakemake_args+=("$1")
            shift
            ;;
    esac
done

if ! command -v snakemake >/dev/null 2>&1; then
    echo "snakemake is not available on PATH" >&2
    exit 127
fi

workflow_seed=$(python - <<PYSEED
from pathlib import Path
import re

config_text = Path("config/config.yaml").read_text()
for line in config_text.splitlines():
    match = re.match(r"^workflow_seed\s*:\s*([^#\s]+)", line)
    if match:
        print(int(match.group(1)))
        break
else:
    raise SystemExit("workflow_seed is missing from config/config.yaml")
PYSEED
)
export SCRNASEQ_PREPROCESS_SEED="${SCRNASEQ_PREPROCESS_SEED:-${workflow_seed}}"
if [[ "${SCRNASEQ_PREPROCESS_SEED}" != "${workflow_seed}" ]]; then
    echo "SCRNASEQ_PREPROCESS_SEED=${SCRNASEQ_PREPROCESS_SEED} does not match config workflow_seed=${workflow_seed}" >&2
    exit 2
fi

profile_args=(--workflow-profile profiles/slurm --profile cannon)
if [[ -n "${TEST_WORKFLOW_PROFILE_ARGS:-}" ]]; then
    # Intended for simple Snakemake profile/executor args.
    read -r -a profile_args <<< "${TEST_WORKFLOW_PROFILE_ARGS}"
fi

common_args=(
    --snakefile workflow/Snakefile
    --configfile config/config.yaml
    --config sampleTable=testdata/samplesheet_test.tsv resultsDir=testdata/results downsampleResultsDir=testdata/results/downsampling
    --use-conda
    "${profile_args[@]}"
)

conda_prefix_args=()
if [[ -n "${snakemake_conda_prefix}" ]]; then
    conda_prefix_args=(--conda-prefix "${snakemake_conda_prefix}")
fi

test_workflow_jobs="${TEST_WORKFLOW_JOBS:-200}"

env_snakemake_args=()
if [[ -n "${TEST_WORKFLOW_SNAKEMAKE_ARGS:-}" ]]; then
    # Intended for simple Snakemake CLI tokens, such as --set-resources entries.
    read -r -a env_snakemake_args <<< "${TEST_WORKFLOW_SNAKEMAKE_ARGS}"
fi

snakemake --unlock "${common_args[@]}" "${conda_prefix_args[@]}"
rm -rf testdata/results

snakemake \
    "${common_args[@]}" \
    "${conda_prefix_args[@]}" \
    --rerun-incomplete \
    --retries 2 \
    --jobs "${test_workflow_jobs}" \
    --latency-wait 600 \
    "${env_snakemake_args[@]}" \
    "${extra_snakemake_args[@]}"
