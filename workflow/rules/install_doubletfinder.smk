localrules: install_doubletfinder

rule install_doubletfinder:
    output:
        f"{RESULTS_DIR}/doubletfinder_installed.txt"
    log:
        f"{RESULTS_DIR}/logs/install_doubletfinder.log"
    conda:
        "../envs/doubletfinder.yml"
    shell:
        """
        Rscript -e "remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)" > {log} 2>&1
        touch {output}
        """
