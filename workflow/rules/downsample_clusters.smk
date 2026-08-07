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
        # All 100 replicates run serially in one job (rm + gc between iterations), so the
        # job's peak memory is a SINGLE replicate's SCTransform footprint. Observed per-
        # replicate peak is ~42 GB and does NOT track input rds size (small datasets such as
        # cteleta are among the heaviest), so a flat baseline beats input-scaling here.
        # Wall-time is ~100x the per-replicate cost (worst observed ~370 s/rep -> ~10.2 h).
        # 64 GB / 15 h cover the worst case on attempt 1; OOM/TIMEOUT restarts are not
        # reliably resubmitted, so baselines must not depend on the retry escalation.
        mem_mb=lambda wildcards, attempt: int(64000 * (2 ** (attempt - 1))),
        runtime=lambda wildcards, attempt: int(900 * (2 ** (attempt - 1)))
    shell:
        """
        SCRNASEQ_DOWNSAMPLE_SEED={params.seed} Rscript {input.script} {input.seurat_object} {output.tsv} {params.n_replicates} {params.downsample_rate} > {log} 2>&1
        """
