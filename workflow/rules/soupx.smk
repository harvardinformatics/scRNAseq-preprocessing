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
        seurat_base =  "results/seurat_filtered/filtered_seurat_" + "{sample}" + ".rds"
    output:
        "results/soupx/seurat_soupx_" + "{sample}" + ".rds"
    conda:
        "../envs/soupx.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1)))
    shell:
        """
        Rscript workflow/scripts/soupx.R  {input.filtered} {input.raw} {input.seurat_base} {output}
        """
