#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH -o output.%j
#SBATCH -e error.%j
#SBATCH --job-name=Run_19
#SBATCH --mem=220gb

#The fermi environment first needs to be sourced:
source /project/ckarwin/astro/chris/my_envs/joint_analysis_main.sh

cd $SLURM_SUBMIT_DIR
python run_3ML_master.py
