rule cellbender:
    input:
        tenx_raw_h5_input
    output:
        base=f"{RESULTS_DIR}/cellbender/cellbender_{{sample}}.h5",
        filtered=f"{RESULTS_DIR}/cellbender/cellbender_{{sample}}_filtered.h5"
    log:
        f"{RESULTS_DIR}/logs/cellbender/cellbender_{{sample}}.log"
    params:
        # per-sample directory ONLY for checkpoints / temp
        workdir="scratch/cellbender/{sample}",
        seed=WORKFLOW_SEED
    container:
        "docker://us.gcr.io/broad-dsde-methods/cellbender:latest"
    resources:
        mem_mb = lambda wildcards, attempt: int(50000 * (2 ** (attempt - 1))),
        slurm_partition = "gpu",
        gres = "gpu:1",
        runtime = 2880
    shell:
        r"""
        set -euo pipefail

        # base_dir is the directory from which Snakemake was launched
        exec > {log} 2>&1
        base_dir=$(pwd)
        export PYTHONHASHSEED={params.seed}

        seed_args=()
        if cellbender remove-background --help 2>&1 | grep -q -- "--seed"; then
            seed_args=(--seed {params.seed})
        fi

        input_h5="{input}"
        case "${{input_h5}}" in
            /*) ;;
            *) input_h5="${{base_dir}}/${{input_h5}}" ;;
        esac

        base_output="{output.base}"
        case "${{base_output}}" in
            /*) ;;
            *) base_output="${{base_dir}}/${{base_output}}" ;;
        esac

        filtered_output="{output.filtered}"
        case "${{filtered_output}}" in
            /*) ;;
            *) filtered_output="${{base_dir}}/${{filtered_output}}" ;;
        esac

        expected_filtered_output="${{base_output%.h5}}_filtered.h5"
        if [ "${{expected_filtered_output}}" != "${{filtered_output}}" ]; then
            echo "CellBender derives filtered output from --output as ${{expected_filtered_output}}, but the rule declares ${{filtered_output}}" >&2
            exit 1
        fi

        # Make sure output and scratch dirs exist on the host
        mkdir -p "$(dirname "${{base_output}}")"
        mkdir -p "{params.workdir}"

        # Run CellBender in the per-sample scratch dir so ckpt.tar.gz is unique.
        cd "{params.workdir}"

        cellbender remove-background             --cuda             "${{seed_args[@]}}"             --input "${{input_h5}}"             --output "${{base_output}}"

        test -s "${{base_output}}"
        test -s "${{filtered_output}}"
        """
