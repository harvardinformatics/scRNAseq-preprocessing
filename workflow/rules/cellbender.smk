rule cellbender:
    input:
        raw=tenx_raw_h5_input,
        script="workflow/scripts/cellbender_adaptive_run.py"
    output:
        base=f"{RESULTS_DIR}/cellbender/cellbender_{{sample}}.h5",
        filtered=f"{RESULTS_DIR}/cellbender/cellbender_{{sample}}_filtered.h5",
        report=f"{RESULTS_DIR}/cellbender/cellbender_{{sample}}_report.html",
        status=f"{RESULTS_DIR}/cellbender/cellbender_{{sample}}_adaptive_status.txt"
    log:
        f"{RESULTS_DIR}/logs/cellbender/cellbender_{{sample}}.log"
    params:
        # per-sample directory ONLY for checkpoints / temp (holds the initial/ and rerun/ runs)
        workdir="scratch/cellbender/{sample}",
        seed=WORKFLOW_SEED,
        learning_rate=CELLBENDER_LEARNING_RATE,
        adaptive="true" if CELLBENDER_ADAPTIVE else "false"
    container:
        "docker://us.gcr.io/broad-dsde-methods/cellbender@sha256:093f2caf1ce4acae4541ea45e52ab7b220ca131ec73b4d1f664b85fe12850bae"
    resources:
        mem_mb = lambda wildcards, attempt: int(50000 * (2 ** (attempt - 1))),
        slurm_partition = "gpu",
        gres = "gpu:1",
        runtime = 2880
    shell:
        # The adaptive wrapper runs CellBender (once, or twice with a halved learning rate when the
        # report's automated assessment recommends it), copies the chosen run's outputs to the
        # declared paths, keeps both runs' reports, and writes the per-sample status file. It
        # resolves relative paths against the launch dir, so we do NOT cd before invoking it.
        r"""
        set -euo pipefail
        exec > {log} 2>&1
        export PYTHONHASHSEED={params.seed}

        mkdir -p "$(dirname "{output.base}")"
        mkdir -p "{params.workdir}"

        python "{input.script}" \
            --input "{input.raw}" \
            --output-base "{output.base}" \
            --output-filtered "{output.filtered}" \
            --report "{output.report}" \
            --status "{output.status}" \
            --sample "{wildcards.sample}" \
            --learning-rate "{params.learning_rate}" \
            --adaptive "{params.adaptive}" \
            --seed "{params.seed}" \
            --workdir "{params.workdir}"

        test -s "{output.base}"
        test -s "{output.filtered}"
        test -s "{output.report}"
        test -s "{output.status}"
        """


rule cellbender_adaptive_summary:
    # Aggregate every sample's adaptive-run status into one run-level TSV, reporting which samples
    # needed a re-run and whether halving the learning rate resolved the learning curve.
    input:
        status=CELLBENDER_STATUS,
        script="workflow/scripts/cellbender_adaptive_summary.py"
    output:
        f"{RESULTS_DIR}/cellbender/cellbender_adaptive_summary.tsv"
    log:
        f"{RESULTS_DIR}/logs/cellbender/cellbender_adaptive_summary.log"
    # Pure stdlib aggregation; reuse the CellBender image (already pulled for the `cellbender`
    # rule) rather than maintain a separate conda env just to run Python.
    container:
        "docker://us.gcr.io/broad-dsde-methods/cellbender@sha256:093f2caf1ce4acae4541ea45e52ab7b220ca131ec73b4d1f664b85fe12850bae"
    shell:
        r"""
        set -euo pipefail
        python "{input.script}" --output "{output}" {input.status} > {log} 2>&1
        """
