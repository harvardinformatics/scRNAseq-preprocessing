N_DOWNSAMPLE_REPLICATES = int(config.get("nDownsampleReplicates", 100))


rule downsample_clusters:
    input:
        seurat_object=lambda wildcards: DOWNSAMPLE_INPUTS_BY_TARGET[wildcards.downsample_target],
        script="workflow/scripts/downsample_clusters.R"
    output:
        tsv=f"{DOWNSAMPLE_RESULTS_DIR}/{{downsample_target}}_clusterdownsampling.tsv"
    log:
        f"{DOWNSAMPLE_RESULTS_DIR}/logs/downsample_clusters/{{downsample_target}}.log"
    params:
        downsample_rate=lambda wildcards: float(config.get("downsampleRate", 0.8)),
        seed=lambda wildcards: int(config.get("workflowSeed", config.get("workflow_seed", 12345))),
        n_replicates=N_DOWNSAMPLE_REPLICATES
    conda:
        "../envs/downsample_clusters.yml"
    wildcard_constraints:
        downsample_target=DOWNSAMPLE_TARGET_REGEX
    resources:
        mem_mb=lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1))),
        runtime=lambda wildcards, attempt: int(270 * (2 ** (attempt - 1)))
    shell:
        """
        SCRNASEQ_DOWNSAMPLE_SEED={params.seed} Rscript {input.script} {input.seurat_object} {output.tsv} {params.n_replicates} {params.downsample_rate} > {log} 2>&1
        """
