def input_function:
    return sampleinfo.loc[sampleinfo["sampleid"] ==wildcards.sample,"tenx_datadir"].values[0]

rule tenx2seuratrds:
    input:
        raw = input_function + "/raw_feature_bc_matrix"
        filtered =  input_function + "/filtered_feature_bc_matrix"
    output:
        "raw_seurat" + "{sample"} + ".rds",
        "filtered_seurat" + "{sample"} + ".rds"
    conda:
        "envs/tenx2seuratrds" 
