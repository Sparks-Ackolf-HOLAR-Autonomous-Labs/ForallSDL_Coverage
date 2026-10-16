#!/bin/bash
# Embedded-SiC particle check (x03_design.json, x = 0.3). Guest partition only:
# grn061 is off-limits by operator directive, the DGX is reserved.
# Submit ONLY after x03_design.json is on the ledger and pushed:
#   sbatch --array=0-9 sbatch_x03.sh
#SBATCH --job-name=mt_sic03
#SBATCH --partition=notchpeak-guest
#SBATCH --account=owner-guest
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --requeue
#SBATCH --output=logs/%x_%A_%a.out
set -euo pipefail
module load openmpi/4.1.6 lammps
cd "$SLURM_SUBMIT_DIR"; mkdir -p results logs
read -r TAG R SHELL SIC SEED <<< "$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" points.txt)"
OUT="results/${TAG}_s${SEED}"
mpirun -np "$SLURM_NTASKS" lmp -in in.sic_particle -log "${OUT}.log" -screen none \
  -var x 0.3 -var R "$R" -var shell "$SHELL" -var sic "$SIC" -var seed "$SEED" \
  -var pot SiCGe_GeC.tersoff -var out "$OUT"
python3 analyze_mp.py "${OUT}.log" "${OUT}.profile" "${OUT}.json"
echo 0 > "${OUT}.exit_status"
