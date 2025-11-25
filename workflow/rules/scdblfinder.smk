rule scdblfinder:
    input:
        "results/{decon_method}/seurat_{decon_method}_{empty_method}_{sample}.rds"
    output:
         "results/scdblfinder/seurat_scdblfinder_{decon_method}_{empty_method}_{sample}.rds"
    conda:
        "../envs/scdblfinder.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1)))
    shell:
        """
        Rscript workflow/scripts/scdblfinder.R  {input} {output}
        """
