#!/bin/bash
set -e  # stop on error

echo "=============================="
echo " Starting case execution"
echo "=============================="

# Load OpenFOAM environment (adjust if needed)
source $WM_PROJECT_DIR/bin/tools/RunFunctions

python3 fg_protocol.py

echo "=============================="
echo " Post-processing"
echo "=============================="

python3 postProc_foamQualityScan_plot.py
