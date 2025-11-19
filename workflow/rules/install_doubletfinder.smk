localrules: install_doubletfinder

rule install_doubletfinder:
  output:
    "results/doubletfinder_installed.txt"
  conda:
    "../envs/doubletfinder.yml"
  shell:
    """
    Rscript -e "remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)"
    touch {output}
    """
