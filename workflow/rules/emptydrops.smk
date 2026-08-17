rule emptydrops:
    input:
        data=tenx_raw_matrix_input,
        script="workflow/scripts/emptydrops.R",
        helper="workflow/scripts/silhouette_utils.R"
    output:
        seurat=f"{RESULTS_DIR}/emptydrops/filtered_seurat_emptydrops_{{sample}}.rds",
        matrixdir=directory(f"{RESULTS_DIR}/emptydrops/{{sample}}_emptydrops_filtered_matrix"),
        nclusters=temp(f"{RESULTS_DIR}/emptydrops/filtered_seurat_emptydrops_{{sample}}_nclusters.txt"),
        cluster_ids=temp(f"{RESULTS_DIR}/emptydrops/filtered_seurat_emptydrops_{{sample}}_cluster_ids.txt")
    log:
        f"{RESULTS_DIR}/logs/emptydrops/emptydrops_{{sample}}.log"
    conda:
        "../envs/emptydrops.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(64000 * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(480* (2 ** (attempt - 1)))
    params:
        seed=WORKFLOW_SEED
    shell:
        """
        SCRNASEQ_PREPROCESS_SEED={params.seed} Rscript {input.script} {input.data} {output.seurat} {output.matrixdir} {output.nclusters} {output.cluster_ids} {wildcards.sample} > {log} 2>&1
        """
