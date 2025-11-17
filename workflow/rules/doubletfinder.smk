def input_function(wildcards):
    return f"results/{wildcards.decon_method}/seurat_{wildcards.decon_method}_{wildcards.sample}.rds"


rule doubletfinder:
    input:
        input_function
    output:
         "results/doubletfinder/sce_doubletfinder_{decon_method}_{sample}.rds"
    conda:
        "../envs/doubletfinder.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(12000 * (2 ** (attempt - 1)))
    shell:
        """
        Rscript workflow/scripts/doubletfinder.R  {input} {output}
        """
