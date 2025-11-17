def input_function(wildcards):
    return f"results/{wildcards.decon_method}/seurat_{wildcards.decon_method}_{wildcards.sample}.rds"


rule scdblfinder:
    input:
        input_function
    output:
         "results/scdblfinder/seurat_scdblfinder_{decon_method}_{sample}.rds"
    conda:
        "../envs/scdblfinder.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1)))
    shell:
        """
        Rscript workflow/scripts/scdblfinder.R  {input} {output}
        """
