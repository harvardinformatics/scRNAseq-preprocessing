def input_function(wildcards):
    tenx_dir = sampleinfo.loc[sampleinfo["sampleid"] == wildcards.sample, "tenx_datadir"].values[0]
    input_path = os.path.join(tenx_dir, "filtered_feature_bc_matrix")
    print(f"Input path for sample {wildcards.sample}: {input_path}")
    return input_path

rule tenx2seuratrds:
    input:
        data=input_function,
        script="workflow/scripts/tenx2seuratrds.R"
    output:
        rds="results/seurat_filtered/filtered_seurat_tenx_" + "{sample}" + ".rds",
        markers="results/seurat_filtered/filtered_seurat_tenx_" + "{sample}" + "_markergenes.csv"
    conda:
        "../envs/tenx2seuratrds.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(480* (2 ** (attempt - 1)))
        
    shell:
        "Rscript {input.script}  {input.data} {output.rds} {output.markers}"
