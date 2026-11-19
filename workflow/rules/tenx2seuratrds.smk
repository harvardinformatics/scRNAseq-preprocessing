rule tenx2seuratrds:
    input:
        data=tenx_filtered_matrix_input,
        mt_genes=sample_mt_genes_input,
        script="workflow/scripts/tenx2seuratrds.R",
        helper="workflow/scripts/silhouette_utils.R",
        regression_helper="workflow/scripts/regression_utils.R"
    output:
        rds=f"{RESULTS_DIR}/seurat_filtered/filtered_seurat_tenx_{{sample}}.rds",
        nclusters=temp(f"{RESULTS_DIR}/seurat_filtered/filtered_seurat_tenx_{{sample}}_nclusters.txt"),
        cluster_ids=temp(f"{RESULTS_DIR}/seurat_filtered/filtered_seurat_tenx_{{sample}}_cluster_ids.txt")
    log:
        f"{RESULTS_DIR}/logs/tenx2seuratrds/tenx2seuratrds_{{sample}}.log"
    conda:
        "../envs/scrnaseq_r.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(48000 * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(480* (2 ** (attempt - 1)))
    params:
        seed=WORKFLOW_SEED,
        vars_to_regress=VARS_TO_REGRESS_ARG,
        mt_genes=sample_mt_genes
    shell:
        "SCRNASEQ_PREPROCESS_SEED={params.seed} Rscript {input.script} {input.data} {output.rds} {output.nclusters} {output.cluster_ids} --vars-to-regress={params.vars_to_regress} --mt-genes={params.mt_genes:q} --sample={wildcards.sample} > {log} 2>&1"
