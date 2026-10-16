#!/bin/bash
#SBATCH --job-name=mt_p3m_bench
#SBATCH --partition=notchpeak-gpu-guest
#SBATCH --account=owner-gpu-guest
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --exclude=notch168,notch294,notch299,notch204
#SBATCH --output=logs/%x_%j.out
cd "$SLURM_SUBMIT_DIR"; mkdir -p logs
P=/uufs/chpc.utah.edu/sys/installdir/r8/mace/0.3.16_cu126/.venv/bin/python
nvidia-smi --query-gpu=name --format=csv,noheader
for id in bench_disp bench_clu; do $P mace_rnemd.py $id --points bench_points.json --out bench/$id.json --ckpt-dir bench_ckpt --bench 300; done
