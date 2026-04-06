#!/bin/bash
set -e  # stop on error

echo "=============================="
echo " Starting batch case execution"
echo "=============================="

# Load OpenFOAM environment
source $WM_PROJECT_DIR/bin/tools/RunFunctions

python3 preProc_eval_porosityAndPermeability.py
python3 preProc_generate_Kfield.py

echo "Running blockMesh..."
blockMesh > log.blockMesh

echo "Running impesFoam3ph..."
impesFoam3ph
