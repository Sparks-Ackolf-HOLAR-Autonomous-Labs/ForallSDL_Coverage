#!/bin/bash
# Guest partition only (grn061 and the DGX not used).
#SBATCH --job-name=mt_sigec_m1c
#SBATCH --partition=notchpeak-guest
#SBATCH --account=owner-guest
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --requeue
#SBATCH --array=0-50
#SBATCH --output=logs/%x_%A_%a.out
cd "$SLURM_SUBMIT_DIR"; mkdir -p logs out ckpt
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
P=/uufs/chpc.utah.edu/sys/installdir/r8/mace/0.3.16_cu126/.venv/bin/python
T=$SLURM_ARRAY_TASK_ID
$P mace_relax_m1c.py "$T" && echo 0 > out/task$T.exit_status || echo $? > out/task$T.exit_status
