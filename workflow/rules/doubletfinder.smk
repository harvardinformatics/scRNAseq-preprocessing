rule doubletfinder:
    input:
        "results/doubletfinder_installed.txt",
        "results/{decon_method}/seurat_{decon_method}_{sample}.rds"
    output:
         "results/doubletfinder/seurat_doubletfinder_{decon_method}_{sample}.rds"
    conda:
        "../envs/doubletfinder.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1)))
    shell:
        """
        Rscript workflow/scripts/doubletfinder.R  {input} {output}
        """
