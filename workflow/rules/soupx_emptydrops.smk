def input_function_raw(wildcards):
    tenx_dir = sampleinfo.loc[sampleinfo["sampleid"] == wildcards.sample, "tenx_datadir"].values[0]
    input_path = os.path.join(tenx_dir, "raw_feature_bc_matrix")
    return input_path


rule soupx_emptydrops:
    input:
        raw=input_function_raw,
        filtered = "results/emptydrops/{sample}_emptydrops_filtered_matrix",
        seurat_base = "results/emptydrops/filtered_seurat_emptydrops_{sample}.rds",
        script="workflow/scripts/soupx.R"
    output:
        rds="results/soupx/seurat_soupx_emptydrops_" + "{sample}" + ".rds",
        nclusters="results/soupx/seurat_soupx_emptydrops_" + "{sample}" + "_nclusters.txt",
        cluster_ids="results/soupx/seurat_soupx_emptydrops_" + "{sample}" + "_cluster_ids.txt"
    conda:
        "../envs/soupx.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(480* (2 ** (attempt - 1)))
    shell:
        """
        Rscript {input.script} {input.filtered} {input.raw} {input.seurat_base} {output.rds} {output.nclusters} {output.cluster_ids}
        """
