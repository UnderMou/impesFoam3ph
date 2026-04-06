#!/bin/bash
set -e  # stop on error

echo "=============================="
echo " Starting case execution"
echo "=============================="

# Load OpenFOAM environment (adjust if needed)
source $WM_PROJECT_DIR/bin/tools/RunFunctions

echo "Running blockMesh..."
blockMesh

echo "------------------------------"
echo " simulation run"
echo "------------------------------"

impesFoam2ph

echo "=============================="
echo " Running post-processing"
echo "=============================="

python3 postProc_checkSolution.py
