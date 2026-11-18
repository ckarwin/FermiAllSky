#!/bin/bash
# Slurm job for the inner-Galaxy 3ML fit. Submit from this folder:
#     sbatch slurm.sh
#SBATCH --time=72:00:00
#SBATCH -o output.%j
#SBATCH -e error.%j
#SBATCH --job-name=Run_19
#SBATCH --mem=220gb

# Activate the environment with Fermipy and the diffuse-analysis branches of
# threeML and astromodels:
source /project/ckarwin/astro/chris/my_envs/joint_analysis_main.sh

cd $SLURM_SUBMIT_DIR
python run_3ML_master.py
