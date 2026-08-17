rule cellbender2seurat:
    input:
        data=f"{RESULTS_DIR}/cellbender/cellbender_{{sample}}_filtered.h5",
        script="workflow/scripts/cellbender2seurat.R",
        helper="workflow/scripts/silhouette_utils.R"
    output:
        rds=f"{RESULTS_DIR}/cellbender_fromraw/seurat_cellbender_fromraw_{{sample}}.rds",
        nclusters=temp(f"{RESULTS_DIR}/cellbender_fromraw/seurat_cellbender_fromraw_{{sample}}_nclusters.txt"),
        cluster_ids=temp(f"{RESULTS_DIR}/cellbender_fromraw/seurat_cellbender_fromraw_{{sample}}_cluster_ids.txt")
    log:
        f"{RESULTS_DIR}/logs/cellbender2seurat/cellbender2seurat_{{sample}}.log"
    conda:
        "../envs/tenx2seuratrds.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(48000 * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(480* (2 ** (attempt - 1)))
    params:
        seed=WORKFLOW_SEED
    shell:
        """
        SCRNASEQ_PREPROCESS_SEED={params.seed} Rscript {input.script} {input.data} {output.rds} {output.nclusters} {output.cluster_ids} > {log} 2>&1
        """
