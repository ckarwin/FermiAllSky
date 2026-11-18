#!/bin/bash
# Slurm job for the all-sky isotropic fit. Submit from this folder:
#     sbatch slurm.sh
# `run` = prepare + fit + plot. After a successful preparation you can switch
# to `fit` and/or `plot` (commented lines below) to refit or replot quickly.
#SBATCH --time=72:00:00
#SBATCH -o output.%j
#SBATCH -e error.%j
#SBATCH --job-name=iso
#SBATCH --mem=220gb

# Activate the Fermipy/Fermitools environment:
source /project/ckarwin/astro/chris/my_envs/fermipy.sh

cd $SLURM_SUBMIT_DIR
python -u fit_isotropic.py run config.yaml
#python fit_isotropic.py fit config.yaml
#python fit_isotropic.py plot config.yaml
