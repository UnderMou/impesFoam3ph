#!/bin/bash
set -e  # stop on error

echo "=============================="
echo " Starting batch case execution"
echo "=============================="

# Load OpenFOAM environment
source $WM_PROJECT_DIR/bin/tools/RunFunctions

# List of case directories
dirs=(CaseA CaseB CaseB_noG)

# Store root directory
ROOT_DIR=$(pwd)

for dir in "${dirs[@]}"; do
    echo "===================================="
    echo " Running case in directory: $dir"
    echo "===================================="

    cd "$ROOT_DIR/$dir" || exit 1

    echo "Running blockMesh..."
    blockMesh > log.blockMesh

    echo "Running impesFoam3ph..."
    impesFoam3ph > log.impesFoam3ph

    echo "Finished case: $dir"
    echo ""

done

# Return to root directory
cd "$ROOT_DIR" || exit 1

echo "=============================="
echo " Running global post-processing"
echo "=============================="

python3 postProc_plotProfiles.py

echo "=============================="
echo " All cases finished successfully"
echo "=============================="