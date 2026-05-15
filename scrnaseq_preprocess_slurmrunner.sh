#!/bin/bash
#SBATCH -p sapphire,shared
#SBATCH -e preprocesstest_%A.err
#SBATCH -o preprocesstest_%A.out
#SBATCH -J preprocesstest
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem 8000
#SBATCH -t 72:00:00

module purge
module load python
conda activate my_snakemake


which snakemake
snakemake --version


# The default location for environments built by snakemake is from within the local project directory at ./.snakemake/conda. However,
# to avoid having to rebuild those environments every time a new instance of the workflow is launched, one can specify a location
# outside of the local project directory. Thus for this script, a single command line argument is provided that is the full
# path to that directory. If you do not with to write to a persistent location, you can either specify ./.snakemake/conda. as the command line argument, or simply remove the --conda-prefix switch below. 

PATH_TO_MY_CONDA_ENVS=$1

# NOTE: the cannon profile is adapted specifically for Harvard's cannon cluster which uses SLURM as
# it's job scheduler. On most Linux (and some macOS) setups, the default location for global profiles is
#  ~/.config/snakemake/, in other words, nested within /.config/snakemake in one's home directory.
# In the case of the the cannon profile, Snakemake looks for it at ~/.config/snakemake/cannon.
# SLURM users will have to adapt this profile, and the setup will have to be adapted if one is
# using SGE or LSF job schedulers. The global profile is where default partitions are specified, and
# default resources for job submissions.


snakemake --unlock --snakefile workflow/Snakefile --configfile config/config.yaml --use-conda --workflow-profile profiles/slurm --profile cannon

snakemake --conda-prefix $PATH_TO_MY_CONDA_ENVS --snakefile workflow/Snakefile --rerun-incomplete --retries 2 --jobs 200 --latency-wait 120 --configfile config/config.yaml --use-conda --workflow-profile profiles/slurm --profile cannon
