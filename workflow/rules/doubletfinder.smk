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
        # DoubletFinder's paramSweep (sct=TRUE) synthesizes ~25% artificial doublets and
        # re-runs SCTransform across a pK sweep, so peak memory is a large multiple of the
        # input object (~43x observed: 2.1 GB rds -> 92 GB). A flat baseline is either
        # wasteful for small samples or fatal for large ones; scale by input size (floor 48 GB).
        mem_mb = lambda wildcards, input, attempt: int(max(48000, 64 * input.size_mb) * (2 ** (attempt - 1))),
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
        # DoubletFinder's paramSweep (sct=TRUE) synthesizes ~25% artificial doublets and
        # re-runs SCTransform across a pK sweep, so peak memory is a large multiple of the
        # input object (~43x observed: 2.1 GB rds -> 92 GB). A flat baseline is either
        # wasteful for small samples or fatal for large ones; scale by input size (floor 48 GB).
        mem_mb = lambda wildcards, input, attempt: int(max(48000, 64 * input.size_mb) * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(480* (2 ** (attempt - 1)))
    params:
        seed=WORKFLOW_SEED
    shell:
        """
        SCRNASEQ_PREPROCESS_SEED={params.seed} Rscript {input.script} {input.data} {output.rds} {output.nclusters} {output.cluster_ids} > {log} 2>&1
        """
