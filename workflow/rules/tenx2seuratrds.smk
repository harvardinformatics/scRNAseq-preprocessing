rule tenx2seuratrds:
    input:
        data=tenx_filtered_matrix_input,
        script="workflow/scripts/tenx2seuratrds.R",
        helper="workflow/scripts/silhouette_utils.R"
    output:
        rds=f"{RESULTS_DIR}/seurat_filtered/filtered_seurat_tenx_{{sample}}.rds",
        nclusters=temp(f"{RESULTS_DIR}/seurat_filtered/filtered_seurat_tenx_{{sample}}_nclusters.txt"),
        cluster_ids=temp(f"{RESULTS_DIR}/seurat_filtered/filtered_seurat_tenx_{{sample}}_cluster_ids.txt")
    log:
        f"{RESULTS_DIR}/logs/tenx2seuratrds/tenx2seuratrds_{{sample}}.log"
    conda:
        "../envs/tenx2seuratrds.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(480* (2 ** (attempt - 1)))
    params:
        seed=WORKFLOW_SEED
    shell:
        "SCRNASEQ_PREPROCESS_SEED={params.seed} Rscript {input.script} {input.data} {output.rds} {output.nclusters} {output.cluster_ids} > {log} 2>&1"
