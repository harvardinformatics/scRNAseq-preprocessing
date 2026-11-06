localrules: install_doubletfinder

rule install_doubletfinder:
    output:
        f"{RESULTS_DIR}/doubletfinder_installed.txt"
    log:
        f"{RESULTS_DIR}/logs/install_doubletfinder.log"
    params:
        # DoubletFinder is not on conda: install one fixed commit (v2.0.6, upstream HEAD on
        # 2026-11-06) so every build of the environment gets the same code.
        repo="chris-mcginnis-ucsf/DoubletFinder@1b244d8f0d54b4b1cb4365639931bbb16f01e1cd"
    conda:
        "../envs/scrnaseq_r.yml"
    shell:
        """
        Rscript -e 'remotes::install_github("{params.repo}", force = TRUE, dependencies = FALSE, upgrade = "never"); stopifnot(requireNamespace("DoubletFinder", quietly = TRUE), requireNamespace("Seurat", quietly = TRUE), requireNamespace("igraph", quietly = TRUE))' > {log} 2>&1
        touch {output}
        """
