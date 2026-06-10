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

common_args=(
    --snakefile workflow/Snakefile
    --configfile config/config.yaml
    --config sampleTable=testdata/samplesheet_test.tsv resultsDir=testdata/results
    --use-conda
    --workflow-profile profiles/slurm
    --profile cannon
)

conda_prefix_args=()
if [[ -n "${snakemake_conda_prefix}" ]]; then
    conda_prefix_args=(--conda-prefix "${snakemake_conda_prefix}")
fi

snakemake --unlock "${common_args[@]}" "${conda_prefix_args[@]}"
rm -rf testdata/results

snakemake \
    "${common_args[@]}" \
    "${conda_prefix_args[@]}" \
    --rerun-incomplete \
    --retries 2 \
    --jobs 200 \
    --latency-wait 600 \
    "${extra_snakemake_args[@]}"
