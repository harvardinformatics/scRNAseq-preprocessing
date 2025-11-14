def input_function(wildcards):
    tenx_dir = sampleinfo.loc[sampleinfo["sampleid"] == wildcards.sample, "tenx_datadir"].values[0]
    input_path = os.path.join(tenx_dir, "raw_feature_bc_matrix")
    return input_path

rule emptydrops:
    input:
        input_function
    output:
         "results/emptydrops/seurat_emptydrops_{sample}.rds"
    conda:
        "../envs/emptydrops.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1)))
    shell:
        """
        Rscript workflow/scripts/emptydrops.R  {input} {output}
        """
