rule posthocfilter:
    input:
        "results/{doublet_method}/seurat_{doubletmethod}_{decon_method}_{empty_method}_{sample}.rds"

    output:
        "results/posthocfilter/seurat_posthocfilt_{doubletmethod}_{decon_method}_{empty_method}_{sample}.rds"
    conda:
        "../envs/posthocfilter.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1)))
    shell:
        """
        Rscript workflow/scripts/posthocfilter.R  {input} {output}
        """
