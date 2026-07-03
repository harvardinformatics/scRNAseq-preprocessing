N_DOWNSAMPLE_REPLICATES = int(config.get("nDownsampleReplicates", 100))


rule downsample_cluster_replicate:
    input:
        seurat_object=lambda wildcards: DOWNSAMPLE_INPUTS_BY_TARGET[wildcards.downsample_target],
        script="workflow/scripts/downsample_cluster_replicate.R"
    output:
        tsv=f"{DOWNSAMPLE_RESULTS_DIR}/replicates/{{downsample_target}}_replicate{{replicate}}.tsv"
    log:
        f"{DOWNSAMPLE_RESULTS_DIR}/logs/downsample_cluster_replicate/{{downsample_target}}_replicate{{replicate}}.log"
    params:
        downsample_rate=lambda wildcards: float(config.get("downsampleRate", 0.8)),
        seed=lambda wildcards: int(config.get("workflowSeed", config.get("workflow_seed", 12345)))
    conda:
        "../envs/downsample_clusters.yml"
    wildcard_constraints:
        downsample_target=DOWNSAMPLE_TARGET_REGEX,
        replicate=r"\d+"
    resources:
        mem_mb=lambda wildcards, attempt: int(24000 * (2 ** (attempt - 1))),
        runtime=lambda wildcards, attempt: int(30 * (2 ** (attempt - 1)))
    shell:
        """
        SCRNASEQ_DOWNSAMPLE_SEED={params.seed} Rscript {input.script} {input.seurat_object} {output.tsv} {wildcards.replicate} {params.downsample_rate} > {log} 2>&1
        """


rule downsample_clusters:
    input:
        replicates=lambda wildcards: [
            f"{DOWNSAMPLE_RESULTS_DIR}/replicates/{wildcards.downsample_target}_replicate{replicate}.tsv"
            for replicate in range(1, N_DOWNSAMPLE_REPLICATES + 1)
        ]
    output:
        tsv=f"{DOWNSAMPLE_RESULTS_DIR}/{{downsample_target}}_clusterdownsampling.tsv"
    log:
        f"{DOWNSAMPLE_RESULTS_DIR}/logs/downsample_clusters/{{downsample_target}}.log"
    wildcard_constraints:
        downsample_target=DOWNSAMPLE_TARGET_REGEX
    resources:
        mem_mb=lambda wildcards, attempt: int(4000 * (2 ** (attempt - 1))),
        runtime=lambda wildcards, attempt: int(60 * (2 ** (attempt - 1)))
    shell:
        """
        {{ head -n1 {input.replicates[0]}; for f in {input.replicates}; do tail -n +2 "$f"; done; }} > {output.tsv} 2> {log}
        """
