rule doubletfinder:
    input:
        install_validation="results/doubletfinder_installed.txt",
        data="results/{decon_method}/seurat_{decon_method}_{empty_method}_{sample}.rds",
        script="workflow/scripts/doubletfinder.R"
    output:
        rds="results/doubletfinder/seurat_doubletfinder_{decon_method}_{empty_method}_{sample}.rds",
        markers="results/doubletfinder/seurat_doubletfinder_{decon_method}_{empty_method}_{sample}_markergenes.csv" 
    conda:
        "../envs/doubletfinder.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(480* (2 ** (attempt - 1)))
    shell:
        """
        Rscript {input.script}  {input.data} {output.rds} {output.markers}
        """
