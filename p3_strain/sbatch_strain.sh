#!/bin/bash
# P3 biaxial-strain SiGe campaign, one array task per line of $POINTS:
#   run_id x T eps dir ncz seed
# Guest partition only (guest first, grn061 off-limits).
# Submit: sbatch --array=0-$(($(wc -l < $POINTS)-1)) --export=ALL,POINTS=file sbatch_strain.sh
#SBATCH --job-name=mt_p3s
#SBATCH --partition=notchpeak-guest
#SBATCH --account=owner-guest
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=24G
#SBATCH --time=4:00:00
#SBATCH --requeue
#SBATCH --output=logs/%x_%A_%a.out
set -euo pipefail
module load openmpi/4.1.6 lammps
POT=/uufs/chpc.utah.edu/sys/installdir/r8/lammps/08022023.3/share/lammps/potentials/SiCGe.tersoff
cd "$SLURM_SUBMIT_DIR"
mkdir -p results logs
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$POINTS")
read -r RID X T EPS DIR NCZ SEED <<< "$LINE"
# requeue after preemption reruns from scratch; each attempt is logged
echo "$(date -u +%FT%TZ) job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} restart=${SLURM_RESTART_COUNT:-0} host=$(hostname)" >> "results/${RID}.attempts"
mpirun -np "$SLURM_NTASKS" lmp -in in.mp_strain \
  -log "results/${RID}.log" -screen none \
  -var x "$X" -var T "$T" -var eps "$EPS" -var dir "$DIR" -var ncz "$NCZ" \
  -var seed "$SEED" -var ncxy 4 -var swapn 50 \
  -var pot "$POT" -var out "results/${RID}"
python3 analyze_mp.py "results/${RID}.log" "results/${RID}.profile" "results/${RID}.json"
