#!/bin/bash
set -e  # stop on error

echo "=============================="
echo " Starting case execution"
echo "=============================="

# Load OpenFOAM environment (adjust if needed)
source $WM_PROJECT_DIR/bin/tools/RunFunctions

echo "Generating K field (K)..."
python3 preProc_generate_Kfield.py

echo "Generating well field (qt)..."
python3 preProc_generate_qtfield.py

echo "Running blockMesh..."
blockMesh

echo "------------------------------"
echo " simulation run"
echo "------------------------------"

impesFoam3ph

echo "------------------------------"
echo " post-Processing"
echo "------------------------------"

python3 postProc_ajustTime.py