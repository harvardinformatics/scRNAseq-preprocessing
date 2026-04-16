def input_function(wildcards):
    tenx_dir = sampleinfo.loc[sampleinfo["sampleid"] == wildcards.sample, "tenx_datadir"].values[0]
    input_path = os.path.join(tenx_dir, "raw_feature_bc_matrix")
    return input_path

rule emptydrops:
    input:
        data=input_function,
        script="workflow/scripts/emptydrops.R"
    output:
        seurat="results/emptydrops/filtered_seurat_emptydrops_{sample}.rds",
        matrixdir=directory("results/emptydrops/{sample}_emptydrops_filtered_matrix"),
        markers="results/emptydrops/filtered_seurat_emptydrops_{sample}_markergenes.csv"
    conda:
        "../envs/emptydrops.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(480* (2 ** (attempt - 1)))
    shell:
        """
        Rscript {input.script}  {input.data} {output.seurat} {output.matrixdir} {output.markers}
        """
