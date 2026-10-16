#!/bin/bash
# Runs the full C-corner gate on the DGX (serial LAMMPS, 1 thread).
set -e
cd "$(dirname "$0")"
export LD_LIBRARY_PATH=$HOME/mat_theory_pilot/mpichlib/usr/lib/aarch64-linux-gnu:$HOME/mat_theory_pilot/mpichlib/lib/aarch64-linux-gnu
export OMP_NUM_THREADS=1
PY=$HOME/mat_theory_pilot/venv/bin/python

$PY gate.py
