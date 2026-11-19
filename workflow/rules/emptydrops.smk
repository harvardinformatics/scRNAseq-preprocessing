rule emptydrops:
    input:
        data=tenx_raw_matrix_input,
        mt_genes=sample_mt_genes_input,
        script="workflow/scripts/emptydrops.R",
        helper="workflow/scripts/silhouette_utils.R",
        regression_helper="workflow/scripts/regression_utils.R"
    output:
        seurat=f"{RESULTS_DIR}/emptydrops/filtered_seurat_emptydrops_{{sample}}.rds",
        matrixdir=directory(f"{RESULTS_DIR}/emptydrops/{{sample}}_emptydrops_filtered_matrix"),
        nclusters=temp(f"{RESULTS_DIR}/emptydrops/filtered_seurat_emptydrops_{{sample}}_nclusters.txt"),
        cluster_ids=temp(f"{RESULTS_DIR}/emptydrops/filtered_seurat_emptydrops_{{sample}}_cluster_ids.txt")
    log:
        f"{RESULTS_DIR}/logs/emptydrops/emptydrops_{{sample}}.log"
    conda:
        "../envs/scrnaseq_r.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(64000 * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(480* (2 ** (attempt - 1)))
    params:
        seed=WORKFLOW_SEED,
        vars_to_regress=VARS_TO_REGRESS_ARG,
        mt_genes=sample_mt_genes
    shell:
        """
        SCRNASEQ_PREPROCESS_SEED={params.seed} Rscript {input.script} {input.data} {output.seurat} {output.matrixdir} {output.nclusters} {output.cluster_ids} {wildcards.sample} --vars-to-regress={params.vars_to_regress} --mt-genes={params.mt_genes:q} > {log} 2>&1
        """
