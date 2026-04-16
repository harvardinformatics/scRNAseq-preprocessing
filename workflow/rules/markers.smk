from pathlib import Path


def read_cluster_ids(cluster_ids_path):
    with open(cluster_ids_path) as handle:
        return [line.strip() for line in handle if line.strip()]


checkpoint marker_manifest:
    input:
        cluster_ids="results/{prefix}_cluster_ids.txt"
    output:
        manifest=directory("results/{prefix}_marker_manifest")
    run:
        manifest_dir = Path(output.manifest)
        manifest_dir.mkdir(parents=True, exist_ok=True)

        for existing_path in manifest_dir.iterdir():
            if existing_path.is_file():
                existing_path.unlink()

        for cluster_id in read_cluster_ids(input.cluster_ids):
            (manifest_dir / f"{cluster_id}.txt").write_text(f"{cluster_id}\n")


def marker_chunk_inputs(wildcards):
    manifest_dir = checkpoints.marker_manifest.get(prefix=wildcards.prefix).output.manifest
    cluster_ids = glob_wildcards(str(Path(manifest_dir) / "{cluster}.txt")).cluster
    return expand(
        "results/{prefix}_markergenes_cluster{cluster}.csv",
        prefix=wildcards.prefix,
        cluster=cluster_ids
    )


rule find_markers:
    input:
        rds="results/{prefix}.rds",
        script="workflow/scripts/find_markers.R"
    output:
        "results/{prefix}_markergenes_cluster{cluster}.csv"
    conda:
        "../envs/tenx2seuratrds.yml"
    resources:
        mem_mb=lambda wildcards, attempt: int(12000 * (2 ** (attempt - 1))),
        runtime=lambda wildcards, attempt: int(240 * (2 ** (attempt - 1)))
    shell:
        "Rscript {input.script} {input.rds} {wildcards.cluster} {output}"


rule combine_markers:
    input:
        markers=marker_chunk_inputs,
        script="workflow/scripts/combine_markers.R"
    output:
        "results/{prefix}_markergenes.csv"
    conda:
        "../envs/tenx2seuratrds.yml"
    resources:
        mem_mb=lambda wildcards, attempt: int(4000 * (2 ** (attempt - 1))),
        runtime=lambda wildcards, attempt: int(60 * (2 ** (attempt - 1)))
    shell:
        "Rscript {input.script} {output} {input.markers}"
