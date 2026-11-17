#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH -o output.%j
#SBATCH -e error.%j
#SBATCH --job-name=iso
#SBATCH --mem=220gb

#The fermi environment first needs to be sourced:
source /project/ckarwin/astro/chris/my_envs/fermipy.sh

cd $SLURM_SUBMIT_DIR
python -u fit_isotropic.py run config.yaml
#python fit_isotropic.py fit config.yaml
#python fit_isotropic.py plot config.yaml
