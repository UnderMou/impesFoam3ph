#!/bin/bash
set -e  # stop on error

echo "=============================="
echo " Starting batch case execution"
echo "=============================="

# Load OpenFOAM environment
source $WM_PROJECT_DIR/bin/tools/RunFunctions

# List of case directories
dirs=(noGravity withGravity)

# Store root directory
ROOT_DIR=$(pwd)

for dir in "${dirs[@]}"; do
    echo "===================================="
    echo " Running case in directory: $dir"
    echo "===================================="

    cd "$ROOT_DIR/$dir" || exit 1

    echo "Running blockMesh..."
    blockMesh > log.blockMesh

    echo "Running impesFoam2ph..."
    impesFoam2ph > log.impesFoam2ph

    echo "Finished case: $dir"
    echo ""

done

# Return to root directory
cd "$ROOT_DIR" || exit 1

echo "=============================="
echo " Post-processing"
echo "=============================="

python3 postProc_analytical.py

echo "=============================="
echo " Finished successfully"
echo "=============================="