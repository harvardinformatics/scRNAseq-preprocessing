checkpoint marker_manifest:
    input:
        cluster_ids=f"{RESULTS_DIR}/{{prefix}}_cluster_ids.txt"
    output:
        manifest=temp(directory(f"{RESULTS_DIR}/{{prefix}}_marker_manifest"))
    log:
        f"{RESULTS_DIR}/logs/markers/{{prefix}}_marker_manifest.log"
    conda:
        "../envs/tenx2seuratrds.yml"
    script:
        "../scripts/write_marker_manifest.py"


rule find_markers:
    input:
        rds=f"{RESULTS_DIR}/{{prefix}}.rds",
        script="workflow/scripts/find_markers.R",
        helper="workflow/scripts/silhouette_utils.R"
    output:
        temp(f"{RESULTS_DIR}/{{prefix}}_markergenes_cluster{{cluster}}.csv")
    log:
        f"{RESULTS_DIR}/logs/markers/{{prefix}}_markergenes_cluster{{cluster}}.log"
    conda:
        "../envs/tenx2seuratrds.yml"
    resources:
        mem_mb=lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1))),
        runtime=lambda wildcards, attempt: int(240 * (2 ** (attempt - 1)))
    params:
        seed=WORKFLOW_SEED
    shell:
        "SCRNASEQ_PREPROCESS_SEED={params.seed} Rscript {input.script} {input.rds} {wildcards.cluster} {output} > {log} 2>&1"


# Temp files are left for Snakemake to delete: a job that removes files the main process
# tracks can crash it ("Unable to obtain modification time ... although it existed
# before") when a checkpoint update re-checks every job's files. The manifest is declared
# as an input so Snakemake removes it once the markers are combined; it never removes a
# checkpoint's temp output that no rule lists as an input.
rule combine_markers:
    input:
        markers=marker_chunk_inputs,
        manifest=f"{RESULTS_DIR}/{{prefix}}_marker_manifest",
        script="workflow/scripts/combine_markers.R"
    output:
        f"{RESULTS_DIR}/{{prefix}}_markergenes.csv"
    log:
        f"{RESULTS_DIR}/logs/markers/{{prefix}}_combine_markers.log"
    conda:
        "../envs/tenx2seuratrds.yml"
    resources:
        mem_mb=lambda wildcards, attempt: int(4000 * (2 ** (attempt - 1))),
        runtime=lambda wildcards, attempt: int(60 * (2 ** (attempt - 1)))
    params:
        seed=WORKFLOW_SEED
    shell:
        "SCRNASEQ_PREPROCESS_SEED={params.seed} Rscript {input.script} {output} {input.markers} > {log} 2>&1"
