def input_function(wildcards):
    return f"results/{wildcards.ambientdecon_method}/{wildcards.sample}.txt"

# Output function: distinguish based on method
def output_function(wildcards):
    return f"results/{wildcards.doublet_method}_{wildcards.sample}.txt"

rule scdblfinder:
    input:
        seurat  =  "results/soupx/filtered_seurat_" + "{sample}" + ".rds"
    output:
        "results/soupx/seurat_scdblfinder_" + "{sample}" + ".rds"
    conda:
        "../envs/soupx.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(slurm_profile_config["mem_mb"]['default-resources']['mem_mb'] * (2 ** (attempt - 1)))
    shell:
        """
        Rscript workflow/scripts/soupx.R  {input.filtered} {input.raw} {input.seurat_base} {output}
        """
