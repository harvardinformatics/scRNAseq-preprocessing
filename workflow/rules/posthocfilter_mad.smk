rule posthocfilter_mad:
    input:
        data="results/{doublet_method}/seurat_{doublet_method}_{decon_method}_{empty_method}_{sample}.rds",
        script="workflow/scripts/posthocfilter_mad.R"

    output:
        rds="results/posthocfilter/seurat_posthocfilt_mad_{doublet_method}_{decon_method}_{empty_method}_{sample}.rds",
        markers="results/posthocfilter/seurat_posthocfilt_mad_{doublet_method}_{decon_method}_{empty_method}_{sample}_markergenes.csv"
    conda:
        "../envs/posthocfilter.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(480* (2 ** (attempt - 1)))
    shell:
        """
        Rscript {input.script}  {input.data} {output.rds} {output.markers}
        """
