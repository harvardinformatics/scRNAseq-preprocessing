def input_function(wildcards):
    tenx_dir = sampleinfo.loc[
        sampleinfo["sampleid"] == wildcards.sample,
        "tenx_datadir"
    ].values[0]
    input_path = os.path.join(tenx_dir, "raw_feature_bc_matrix.h5")
    return input_path


rule cellbender:
    input:
        input_function
    output:
        base     = "results/cellbender/cellbender_{sample}.h5",
        filtered = "results/cellbender/cellbender_{sample}_filtered.h5"
    params:
        # per-sample directory ONLY for checkpoints / temp
        workdir = "scratch/cellbender/{sample}"
    container:
        "docker://us.gcr.io/broad-dsde-methods/cellbender:latest"
    resources:
        mem_mb = lambda wildcards, attempt: int(50000 * (2 ** (attempt - 1))),
        slurm_partition = "gpu",
        gres = "gpu:1",
        runtime = 2880
    shell:
        r"""
        # base_dir is the directory from which Snakemake was launched
        base_dir=$(pwd)

        # Make sure output and scratch dirs exist on the host
        mkdir -p "$(dirname "{output.base}")"
        mkdir -p "{params.workdir}"

        # Run cellbender in the per-sample scratch dir so ckpt.tar.gz is unique
        cd "{params.workdir}"

        cellbender remove-background \
            --cuda \
            --input "{input}" \
            --output "${{base_dir}}/{output.base}"

        # If you actually generate a separate filtered file, add that command here,
        # and ensure it writes to ${{base_dir}}/{output.filtered}
        """
