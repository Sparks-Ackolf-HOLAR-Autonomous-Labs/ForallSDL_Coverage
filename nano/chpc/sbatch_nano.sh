#!/bin/bash
# Guest partition only (standing directive: guest/free first, grn061 off-limits).
#SBATCH --job-name=mt_nano
#SBATCH --partition=notchpeak-guest
#SBATCH --account=owner-guest
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=4:00:00
#SBATCH --requeue
#SBATCH --array=0-8
#SBATCH --output=logs/%x_%A_%a.out
XS=(0 0.1 0.15 0.2 0.3 0.5 0.7 0.85 1)
cd "$SLURM_SUBMIT_DIR"; mkdir -p logs
bash run_nano_x.sh "${XS[$SLURM_ARRAY_TASK_ID]}"
