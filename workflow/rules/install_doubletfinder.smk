localrules: install_doubletfinder

rule install_doubletfinder:
  input:
    # whatever required
  output:
    touch("envs/doubletfinder_installed.txt")
  conda:
    "envs/doubletfinder.yaml"
  shell:
    """
    Rscript -e "remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)"
    touch {output}
    """
