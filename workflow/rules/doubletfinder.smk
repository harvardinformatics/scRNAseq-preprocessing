def input_function(wildcards):
    return f"results/{wildcards.decon_method}/seurat_{wildcards.decon_method}_{wildcards.sample}.rds"


rule doubletfinder:
    input:
        input_function
    output:
         "results/doubletfinder/seurat_doubletfinder_{decon_method}_{sample}.rds"
    conda:
        "../envs/doubletfinder.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1)))
    shell:
        """
        which R
        which Rscript
        echo $PATH
        Rscript workflow/scripts/doubletfinder.R  {input} {output}
        """
