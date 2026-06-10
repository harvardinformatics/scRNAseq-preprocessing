rule doubletfinder:
    input:
        install_validation=f"{RESULTS_DIR}/doubletfinder_installed.txt",
        data=f"{RESULTS_DIR}/{{decon_method}}/seurat_{{decon_method}}_{{empty_method}}_{{sample}}.rds",
        script="workflow/scripts/doubletfinder.R",
        helper="workflow/scripts/silhouette_utils.R"
    output:
        rds=f"{RESULTS_DIR}/doubletfinder/seurat_doubletfinder_{{decon_method}}_{{empty_method}}_{{sample}}.rds",
        nclusters=temp(f"{RESULTS_DIR}/doubletfinder/seurat_doubletfinder_{{decon_method}}_{{empty_method}}_{{sample}}_nclusters.txt"),
        cluster_ids=temp(f"{RESULTS_DIR}/doubletfinder/seurat_doubletfinder_{{decon_method}}_{{empty_method}}_{{sample}}_cluster_ids.txt")
    log:
        f"{RESULTS_DIR}/logs/doubletfinder/doubletfinder_{{decon_method}}_{{empty_method}}_{{sample}}.log"
    conda:
        "../envs/doubletfinder.yml"
    wildcard_constraints:
        decon_method="soupx",
        empty_method="tenx|emptydrops"
    resources:
        mem_mb = lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(480* (2 ** (attempt - 1)))
    params:
        seed=WORKFLOW_SEED
    shell:
        """
        SCRNASEQ_PREPROCESS_SEED={params.seed} Rscript {input.script} {input.data} {output.rds} {output.nclusters} {output.cluster_ids} > {log} 2>&1
        """


rule doubletfinder_cellbender:
    input:
        install_validation=f"{RESULTS_DIR}/doubletfinder_installed.txt",
        data=f"{RESULTS_DIR}/cellbender_fromraw/seurat_cellbender_fromraw_{{sample}}.rds",
        script="workflow/scripts/doubletfinder.R",
        helper="workflow/scripts/silhouette_utils.R"
    output:
        rds=f"{RESULTS_DIR}/doubletfinder/seurat_doubletfinder_cellbender_fromraw_{{sample}}.rds",
        nclusters=temp(f"{RESULTS_DIR}/doubletfinder/seurat_doubletfinder_cellbender_fromraw_{{sample}}_nclusters.txt"),
        cluster_ids=temp(f"{RESULTS_DIR}/doubletfinder/seurat_doubletfinder_cellbender_fromraw_{{sample}}_cluster_ids.txt")
    log:
        f"{RESULTS_DIR}/logs/doubletfinder/doubletfinder_cellbender_fromraw_{{sample}}.log"
    conda:
        "../envs/doubletfinder.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(480* (2 ** (attempt - 1)))
    params:
        seed=WORKFLOW_SEED
    shell:
        """
        SCRNASEQ_PREPROCESS_SEED={params.seed} Rscript {input.script} {input.data} {output.rds} {output.nclusters} {output.cluster_ids} > {log} 2>&1
        """
