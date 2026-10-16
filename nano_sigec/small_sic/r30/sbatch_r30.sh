#!/bin/bash
# Fixed-volume-fraction SiC radius check (r30_design.json, x = 0.2).
# Submit ONLY after r30_design.json is on the ledger and pushed:
#   sbatch --array=0-7 sbatch_r30.sh
#SBATCH --job-name=mt_sicr30
#SBATCH --partition=notchpeak-guest
#SBATCH --account=owner-guest
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --requeue
#SBATCH --output=logs/%x_%A_%a.out
set -euo pipefail
module load openmpi/4.1.6 lammps
cd "$SLURM_SUBMIT_DIR"; mkdir -p results logs
read -r TAG R SHELL SIC NCXY NCZ SEED <<< "$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" points.txt)"
OUT="results/${TAG}_s${SEED}"
echo "$(date -u +%FT%TZ) $(hostname) ${SLURM_JOB_ID}" >> "${OUT}.attempts"
mpirun -np "$SLURM_NTASKS" lmp -in in.sic_particle -log "${OUT}.log" -screen none \
  -var x 0.2 -var R "$R" -var shell "$SHELL" -var sic "$SIC" -var seed "$SEED" \
  -var ncxy "$NCXY" -var ncz "$NCZ" -var pot SiCGe_GeC.tersoff -var out "$OUT"
python3 analyze_mp.py "${OUT}.log" "${OUT}.profile" "${OUT}.json"
echo 0 > "${OUT}.exit_status"
