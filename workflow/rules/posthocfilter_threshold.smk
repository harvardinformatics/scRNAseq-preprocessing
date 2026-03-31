rule posthocfilter_threshold:
    input:
        "results/{doublet_method}/seurat_{doublet_method}_{decon_method}_{empty_method}_{sample}.rds"

    output:
        rds="results/posthocfilter/seurat_posthocfilt_threshold_{doublet_method}_{decon_method}_{empty_method}_{sample}.rds",
        markers="results/posthocfilter/seurat_posthocfilt_threshold_{doublet_method}_{decon_method}_{empty_method}_{sample}_markergenes.csv"
    conda:
        "../envs/posthocfilter.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(480* (2 ** (attempt - 1)))
    params:
        min_numfeatures = config["min_nfeature"],
        min_umicount = config["min_ncount"],
        max_mtdna_pcent = config["max_mtdna"]
    shell:
        """
        Rscript workflow/scripts/posthocfilter_threshold.R  {input} {output.rds} \
        {params.min_numfeatures} {params.min_umicount} {params.max_mtdna_pcent} {output.markers} 
        """
