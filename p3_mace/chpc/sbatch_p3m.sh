#!/bin/bash
# P3 hidden-subset study (../p3_mace_design.json, ledger seq 64-65). Guest GPUs only
# (grn061 off-limits; DGX reserved). One task per line of ids.txt. Checkpointed.
#   sbatch --array=0-88 --gres=gpu:1 sbatch_p3m.sh
#SBATCH --job-name=mt_p3m
#SBATCH --partition=notchpeak-gpu-guest
#SBATCH --account=owner-gpu-guest
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --requeue
#SBATCH --exclude=notch168,notch294,notch299,notch204
#SBATCH --output=logs/%x_%A_%a.out
cd "$SLURM_SUBMIT_DIR"; mkdir -p logs out ckpt
ID=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" ids.txt)
[ -s out/$ID.json ] && exit 0
nvidia-smi --query-gpu=name --format=csv,noheader || exit 3
P=/uufs/chpc.utah.edu/sys/installdir/r8/mace/0.3.16_cu126/.venv/bin/python
$P mace_rnemd.py "$ID" --points points.json --out out/$ID.json --ckpt-dir ckpt
