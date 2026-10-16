#!/bin/bash
#SBATCH --job-name=aif_sr
#SBATCH --partition=notchpeak-guest
#SBATCH --account=owner-guest
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=2:00:00
#SBATCH --requeue
#SBATCH -o aif_%j.out
module load python/3.11.7
cd $HOME/mat_theory/aif
$HOME/mat_theory/aif_venv/bin/python aif_run.py
