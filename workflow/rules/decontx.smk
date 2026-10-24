# DecontX ambient-RNA decontamination (campbio/decontX), the R analog of the soupx arm. decontX is
# run with the raw matrix as background (the SoupX-equivalent mode): x = the called-cell (filtered)
# matrix, background = the raw matrix, from which decontX fixes the ambient profile to the
# empty-droplet gene totals (it drops the called barcodes from the background automatically). The
# arm fans over the empty_method (tenx | emptydrops) cell set and uses the same 2-token
# {decon_method}_{empty_method} output naming as soupx, so the generic doubletfinder / scdblfinder /
# posthocfilter rules carry it downstream once their decon_method wildcard_constraints include
# decontx_raw. The filtered matrix and the base Seurat object (whose clusters decontX reuses as its
# cell populations, z) are selected by empty_method via the common.smk helpers.

rule decontx_raw:
    input:
        filtered=decontx_cells_matrix_input,
        raw=tenx_raw_matrix_input,
        seurat_base=decontx_seurat_base_input,
        script="workflow/scripts/decontx.R",
        helper="workflow/scripts/silhouette_utils.R"
    output:
        rds=f"{RESULTS_DIR}/decontx_raw/seurat_decontx_raw_{{empty_method}}_{{sample}}.rds",
        nclusters=temp(f"{RESULTS_DIR}/decontx_raw/seurat_decontx_raw_{{empty_method}}_{{sample}}_nclusters.txt"),
        cluster_ids=temp(f"{RESULTS_DIR}/decontx_raw/seurat_decontx_raw_{{empty_method}}_{{sample}}_cluster_ids.txt")
    log:
        f"{RESULTS_DIR}/logs/decontx/decontx_raw_{{empty_method}}_{{sample}}.log"
    conda:
        "../envs/decontx.yml"
    wildcard_constraints:
        empty_method="tenx|emptydrops"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: int(64000 * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(480* (2 ** (attempt - 1)))
    params:
        seed=WORKFLOW_SEED
    shell:
        """
        SCRNASEQ_PREPROCESS_SEED={params.seed} DECONTX_THREADS={threads} Rscript {input.script} {input.filtered} {input.raw} {input.seurat_base} {output.rds} {output.nclusters} {output.cluster_ids} > {log} 2>&1
        """
