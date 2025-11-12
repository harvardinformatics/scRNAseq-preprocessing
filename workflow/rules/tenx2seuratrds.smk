def input_function(wildcards):
    tenx_dir = sampleinfo.loc[sampleinfo["sampleid"] == wildcards.sample, "tenx_datadir"].values[0]
    input_path = os.path.join(tenx_dir, "filtered_feature_bc_matrix")
    print(f"Input path for sample {wildcards.sample}: {input_path}")
    return input_path

rule tenx2seuratrds:
    input:
        input_function
    output:
        "results/seurat_filtered/filtered_seurat_" + "{sample}" + ".rds"
    conda:
        "../envs/tenx2seuratrds.yml" 
    shell:
        "Rscript workflow/scripts/tenx2seuratrds.R  {input} {output}"
