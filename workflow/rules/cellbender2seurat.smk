rule cellbender2seurat:
    input:
        "results/cellbender/cellbender_{sample}_filtered.h5"
    output:
         "results/cellbender/seurat_cellbender_fromraw_{sample}.rds"
    conda:
        "../envs/tenx2seuratrds.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1)))
    shell:
        """
        Rscript workflow/scripts/cellbender2seurat.R {input} {output}
        """
