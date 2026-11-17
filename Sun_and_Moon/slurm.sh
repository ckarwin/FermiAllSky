#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH -o output.%j
#SBATCH -e error.%j
#SBATCH --job-name=Run_19
#SBATCH --mem=220gb

#The fermi environment first needs to be sourced:
source /project/ckarwin/astro/chris/my_envs/joint_analysis_main.sh

cd $SLURM_SUBMIT_DIR
python -u generate_sun_moon_templates.py /project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Isotropic/isotropic_fit/config.yaml --events /project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Isotropic/isotropic_fit/output/fermipy/ft1_00.fits --ltcube /project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Isotropic/isotropic_fit/output/fermipy/ltcube_00.fits --outdir solar_parallel --finalize-array
