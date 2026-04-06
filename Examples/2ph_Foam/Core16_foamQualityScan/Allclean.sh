#!/bin/bash

echo "=============================="
echo " Cleaning case directories"
echo "=============================="

# Remove files (ignore if they don't exist)
rm -f DP_openFOAM.pdf steadyState_muApp_fg_openFOAM.pdf

# Remove directories
rm -rf fg60.0 fg70.0 fg75.0 fg80.0 fg90.0

# Clean specific case using OpenFOAM utility
if [ -d "fg50.0" ]; then
    echo "Cleaning fg50.0 with foamList..."
    foamListTimes -rm -case fg50.0
fi

echo "=============================="
echo " Cleanup finished"
echo "=============================="