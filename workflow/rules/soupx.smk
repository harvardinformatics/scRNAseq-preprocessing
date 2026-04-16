def input_function_filtered(wildcards):
    tenx_dir = sampleinfo.loc[sampleinfo["sampleid"] == wildcards.sample, "tenx_datadir"].values[0]
    input_path = os.path.join(tenx_dir, "filtered_feature_bc_matrix")
    return input_path

def input_function_raw(wildcards):
    tenx_dir = sampleinfo.loc[sampleinfo["sampleid"] == wildcards.sample, "tenx_datadir"].values[0]
    input_path = os.path.join(tenx_dir, "raw_feature_bc_matrix")
    return input_path


rule soupx:
    input:
        raw=input_function_raw,
        filtered =  input_function_filtered,
        seurat_base =  "results/seurat_filtered/filtered_seurat_tenx_" + "{sample}" + ".rds",
        script="workflow/scripts/soupx.R"
    output:
        rds="results/soupx/seurat_soupx_tenx_" + "{sample}" + ".rds",
        markers="results/soupx/seurat_soupx_tenx_" + "{sample}" + "_markergenes.csv"
    conda:
        "../envs/soupx.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(480* (2 ** (attempt - 1)))
    shell:
        """
        Rscript {input.script}  {input.filtered} {input.raw} {input.seurat_base} {output.rds} {output.markers}
        """
