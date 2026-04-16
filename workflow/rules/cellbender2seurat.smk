rule cellbender2seurat:
    input:
        data="results/cellbender/cellbender_{sample}_filtered.h5",
        script="workflow/scripts/cellbender2seurat.R"
    output:
         rds="results/cellbender/seurat_cellbender_fromraw_{sample}.rds",
         nclusters="results/cellbender/seurat_cellbender_fromraw_{sample}_nclusters.txt",
         cluster_ids="results/cellbender/seurat_cellbender_fromraw_{sample}_cluster_ids.txt"
    conda:
        "../envs/tenx2seuratrds.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(480* (2 ** (attempt - 1)))
    shell:
        """
        Rscript {input.script} {input.data} {output.rds} {output.nclusters} {output.cluster_ids}
        """
