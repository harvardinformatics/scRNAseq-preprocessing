rule soupx:
    input:
        raw=tenx_raw_matrix_input,
        filtered=tenx_filtered_matrix_input,
        seurat_base=f"{RESULTS_DIR}/seurat_filtered/filtered_seurat_tenx_{{sample}}.rds",
        script="workflow/scripts/soupx.R",
        helper="workflow/scripts/silhouette_utils.R"
    output:
        rds=f"{RESULTS_DIR}/soupx/seurat_soupx_tenx_{{sample}}.rds",
        nclusters=temp(f"{RESULTS_DIR}/soupx/seurat_soupx_tenx_{{sample}}_nclusters.txt"),
        cluster_ids=temp(f"{RESULTS_DIR}/soupx/seurat_soupx_tenx_{{sample}}_cluster_ids.txt")
    log:
        f"{RESULTS_DIR}/logs/soupx/soupx_tenx_{{sample}}.log"
    conda:
        "../envs/soupx.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(480* (2 ** (attempt - 1)))
    params:
        seed=WORKFLOW_SEED
    shell:
        """
        SCRNASEQ_PREPROCESS_SEED={params.seed} Rscript {input.script} {input.filtered} {input.raw} {input.seurat_base} {output.rds} {output.nclusters} {output.cluster_ids} > {log} 2>&1
        """
