rule scdblfinder:
    input:
        data="results/{decon_method}/seurat_{decon_method}_{empty_method}_{sample}.rds",
        script="workflow/scripts/scdblfinder.R"
    output:
        rds="results/scdblfinder/seurat_scdblfinder_{decon_method}_{empty_method}_{sample}.rds",
        nclusters="results/scdblfinder/seurat_scdblfinder_{decon_method}_{empty_method}_{sample}_nclusters.txt",
        cluster_ids="results/scdblfinder/seurat_scdblfinder_{decon_method}_{empty_method}_{sample}_cluster_ids.txt"
    conda:
        "../envs/scdblfinder.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(480* (2 ** (attempt - 1)))
    shell:
        """
        Rscript {input.script} {input.data} {output.rds} {output.nclusters} {output.cluster_ids}
        """
