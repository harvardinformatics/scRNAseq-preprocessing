# CellBender's cell calls without its background removal (config: cellbender_cell_calls). The raw
# counts of the barcodes in CellBender's *_filtered.h5 are built into a Seurat object by the
# CellRanger rule's own script, so CellBender's cell calling compares directly with CellRanger's
# (seurat_filtered) and emptyDrops' (emptydrops) before any ambient-RNA decontamination.
rule cellbendercells_matrix:
    input:
        raw=tenx_raw_matrix_input,
        cellbender=f"{RESULTS_DIR}/cellbender/cellbender_{{sample}}_filtered.h5",
        script="workflow/scripts/cellbendercells_matrix.R"
    output:
        matrixdir=directory(f"{RESULTS_DIR}/cellbendercells/{{sample}}_cellbendercells_filtered_matrix")
    log:
        f"{RESULTS_DIR}/logs/cellbendercells_matrix/cellbendercells_matrix_{{sample}}.log"
    conda:
        "../envs/emptydrops.yml"
    resources:
        # Reads the whole raw matrix, as emptydrops does.
        mem_mb = lambda wildcards, attempt: int(64000 * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(240 * (2 ** (attempt - 1)))
    shell:
        "Rscript {input.script} {input.raw} {input.cellbender} {output.matrixdir} > {log} 2>&1"


rule cellbendercells2seurat:
    input:
        data=f"{RESULTS_DIR}/cellbendercells/{{sample}}_cellbendercells_filtered_matrix",
        script="workflow/scripts/tenx2seuratrds.R",
        helper="workflow/scripts/silhouette_utils.R"
    output:
        rds=f"{RESULTS_DIR}/cellbendercells/filtered_seurat_cellbendercells_{{sample}}.rds",
        nclusters=temp(f"{RESULTS_DIR}/cellbendercells/filtered_seurat_cellbendercells_{{sample}}_nclusters.txt"),
        cluster_ids=temp(f"{RESULTS_DIR}/cellbendercells/filtered_seurat_cellbendercells_{{sample}}_cluster_ids.txt")
    log:
        f"{RESULTS_DIR}/logs/cellbendercells2seurat/cellbendercells2seurat_{{sample}}.log"
    conda:
        "../envs/tenx2seuratrds.yml"
    resources:
        mem_mb = lambda wildcards, attempt: int(48000 * (2 ** (attempt - 1))),
        runtime = lambda wildcards, attempt: int(480* (2 ** (attempt - 1)))
    params:
        seed=WORKFLOW_SEED
    shell:
        "SCRNASEQ_PREPROCESS_SEED={params.seed} Rscript {input.script} {input.data} {output.rds} {output.nclusters} {output.cluster_ids} > {log} 2>&1"
