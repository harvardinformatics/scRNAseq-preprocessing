def input_function(wildcards):
    return f"results/{wildcards.decon_method}/seurat_{wildcards.decon_method}_{wildcards.sample}.rds"


rule scdblfinder:
    input:
        input_function
    output:
         "results/{doublet_method}/seurat_{doublet_method}_{decon_method}_{sample}.rds"
    conda:
        "../envs/scdblfinder.yml"
    shell:
        """
        Rscript workflow/scripts/scdblfinder.R  {input} {output}
        """
