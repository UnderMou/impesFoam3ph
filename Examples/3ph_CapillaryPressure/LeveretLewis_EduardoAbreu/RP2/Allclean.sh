#!/bin/bash

echo "=============================="
echo " Cleaning case directories"
echo "=============================="

# Remove files (ignore if they don't exist)
rm -f DP_openFOAM.pdf steadyState_muApp_fg_openFOAM.pdf

# Remove directories
dirs=(eps_0.001  eps_0.005  eps_0.01  eps_0.02)

# Store root directory
ROOT_DIR=$(pwd)

for dir in "${dirs[@]}"; do

    cd "$ROOT_DIR/$dir" || exit 1

    foamListTimes -rm

done

# Return to root directory
cd "$ROOT_DIR" || exit 1

echo "=============================="
echo " Cleanup finished"
echo "=============================="