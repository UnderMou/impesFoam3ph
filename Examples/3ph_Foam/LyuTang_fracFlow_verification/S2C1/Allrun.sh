#!/bin/bash
set -e  # stop on error

echo "=============================="
echo " Starting batch case execution"
echo "=============================="

# Load OpenFOAM environment
source $WM_PROJECT_DIR/bin/tools/RunFunctions

echo "Running blockMesh..."
blockMesh > log.blockMesh

echo "Running impesFoam3ph..."
impesFoam3ph

echo "=============================="
echo " Post-processing"
echo "=============================="

python3 postProc_plotProfiles.py

echo "=============================="
echo " Finished successfully"
echo "=============================="