rule posthocfilter_mad:
    input:
        data=f"{RESULTS_DIR}/{{doublet_method}}/seurat_{{doublet_method}}_{{decon_method}}_{{empty_method}}_{{sample}}.rds",
        script="workflow/scripts/posthocfilter_mad.R",
        helper="workflow/scripts/silhouette_utils.R"
    output:
        rds=f"{RESULTS_DIR}/posthocfilter/seurat_posthocfilt_mad_{{doublet_method}}_{{decon_method}}_{{empty_method}}_{{sample}}.rds",
        nclusters=temp(f"{RESULTS_DIR}/posthocfilter/seurat_posthocfilt_mad_{{doublet_method}}_{{decon_method}}_{{empty_method}}_{{sample}}_nclusters.txt"),
        cluster_ids=temp(f"{RESULTS_DIR}/posthocfilter/seurat_posthocfilt_mad_{{doublet_method}}_{{decon_method}}_{{empty_method}}_{{sample}}_cluster_ids.txt")
    log:
        f"{RESULTS_DIR}/logs/posthocfilter/posthocfilter_mad_{{doublet_method}}_{{decon_method}}_{{empty_method}}_{{sample}}.log"
    conda:
        "../envs/posthocfilter.yml"
    wildcard_constraints:
        doublet_method="doubletfinder|scdblfinder",
        decon_method="soupx",
        empty_method="tenx|emptydrops"
    resources:
        mem_mb = lambda wildcards, attempt: int(48000 * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(480* (2 ** (attempt - 1)))
    params:
        seed=WORKFLOW_SEED
    shell:
        """
        SCRNASEQ_PREPROCESS_SEED={params.seed} Rscript {input.script} {input.data} {output.rds} {output.nclusters} {output.cluster_ids} > {log} 2>&1
        """


rule posthocfilter_mad_cellbender:
    input:
        data=f"{RESULTS_DIR}/{{doublet_method}}/seurat_{{doublet_method}}_cellbender_fromraw_{{sample}}.rds",
        script="workflow/scripts/posthocfilter_mad.R",
        helper="workflow/scripts/silhouette_utils.R"
    output:
        rds=f"{RESULTS_DIR}/posthocfilter/seurat_posthocfilt_mad_{{doublet_method}}_cellbender_fromraw_{{sample}}.rds",
        nclusters=temp(f"{RESULTS_DIR}/posthocfilter/seurat_posthocfilt_mad_{{doublet_method}}_cellbender_fromraw_{{sample}}_nclusters.txt"),
        cluster_ids=temp(f"{RESULTS_DIR}/posthocfilter/seurat_posthocfilt_mad_{{doublet_method}}_cellbender_fromraw_{{sample}}_cluster_ids.txt")
    log:
        f"{RESULTS_DIR}/logs/posthocfilter/posthocfilter_mad_{{doublet_method}}_cellbender_fromraw_{{sample}}.log"
    conda:
        "../envs/posthocfilter.yml"
    wildcard_constraints:
        doublet_method="doubletfinder|scdblfinder"
    resources:
        mem_mb = lambda wildcards, attempt: int(48000 * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(480* (2 ** (attempt - 1)))
    params:
        seed=WORKFLOW_SEED
    shell:
        """
        SCRNASEQ_PREPROCESS_SEED={params.seed} Rscript {input.script} {input.data} {output.rds} {output.nclusters} {output.cluster_ids} > {log} 2>&1
        """
