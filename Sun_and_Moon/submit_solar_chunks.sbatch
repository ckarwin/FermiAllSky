#!/bin/bash
# Slurm array job: each task computes one Sun or Moon livetime chunk.
# Usage (after generate_sun_moon_templates.py --prepare-array; see README.md):
#     mkdir -p palmetto_logs
#     sbatch --array="$(cat solar_parallel/array_range.txt)" submit_solar_chunks.sbatch solar_parallel/array_plan.json
#SBATCH --time=04:00:00
#SBATCH --output=palmetto_logs/solar.%A_%a.out
#SBATCH --error=palmetto_logs/solar.%A_%a.err
#SBATCH --job-name=solar_chunks
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
# Set the array range at submission: sbatch --array=0-LAST%1000 ...
# Add your account/partition directives here if required.
set -euo pipefail
# Activate your Fermitools environment BEFORE submitting this script.
# 64G per task is a conservative example, not a measured requirement.
# Adjust memory after checking your jobs; do not assume 3G is sufficient.
sleep "$((RANDOM % 60))"
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg
"${PYTHON_BIN:-python}" -u solar_array.py chunk "${1:?Pass the path to array_plan.json}" "$SLURM_ARRAY_TASK_ID"
