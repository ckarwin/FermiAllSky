#!/bin/bash
# Final step of the Sun/Moon template production (see README.md):
# sums the completed gtltcubesun array chunks and writes sun_template.fits and
# moon_template.fits. Submit only after every array task has finished:
#     sbatch slurm.sh
#SBATCH --time=72:00:00
#SBATCH -o output.%j
#SBATCH -e error.%j
#SBATCH --job-name=Run_19
#SBATCH --mem=220gb

# Activate the Fermitools environment:
source /project/ckarwin/astro/chris/my_envs/joint_analysis_main.sh

cd $SLURM_SUBMIT_DIR
python -u generate_sun_moon_templates.py /project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Isotropic/isotropic_fit/config.yaml --events /project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Isotropic/isotropic_fit/output/fermipy/ft1_00.fits --ltcube /project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Isotropic/isotropic_fit/output/fermipy/ltcube_00.fits --outdir solar_parallel --finalize-array
