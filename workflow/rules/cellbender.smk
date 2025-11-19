def input_function(wildcards):
    tenx_dir = sampleinfo.loc[sampleinfo["sampleid"] == wildcards.sample, "tenx_datadir"].values[0]
    input_path = os.path.join(tenx_dir, "raw_feature_bc_matrix.h5")
    return input_path

rule cellbender:
    input:
        input_function
    output:
         "results/cellbender/seurat_cellbender_{sample}.rds"
    conda:
        "../envs/cellbender.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(50000 * (2 ** (attempt - 1))),
        slurm_partition="gpu",    
        gres="gpu:1",
        runtime=2880   
    shell:
        """
        cellbender remove-background --cuda --input {input} --output results/cellbender/cellbender_{wildcards.sample}.h5
        Rscript workflow/script/cellbender.R results/cellbender/cellbender_{wildcards.sample}.h5 {output}
        """
