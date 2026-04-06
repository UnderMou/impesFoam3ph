#!/bin/bash
set -e  # stop on error

echo "=============================="
echo " Starting case execution"
echo "=============================="

# Load OpenFOAM environment (adjust if needed)
source $WM_PROJECT_DIR/bin/tools/RunFunctions

cp system/controlDict.orig system/controlDict

echo "Running blockMesh..."
blockMesh

echo "------------------------------"
echo " First simulation run"
echo "------------------------------"

impesFoam2ph

echo "------------------------------"
echo " Modifying controlDict"
echo "------------------------------"

CONTROL_DICT="system/controlDict"

# --- New values ---
newDeltaT=1e1
newEndTime=2e6
newWriteInterval=1e3

# Replace entries using sed
sed -i "s/^deltaT.*/deltaT         ${newDeltaT};/" $CONTROL_DICT
sed -i "s/^endTime.*/endTime       ${newEndTime};/" $CONTROL_DICT
sed -i "s/^writeInterval.*/writeInterval ${newWriteInterval};/" $CONTROL_DICT

echo "Updated controlDict:"
grep -E "deltaT|endTime|writeInterval" $CONTROL_DICT

echo "------------------------------"
echo " Second simulation run"
echo "------------------------------"

impesFoam2ph

echo "------------------------------"
echo " Running post-processing"
echo "------------------------------"

python3 postProc_ajustTime.py
python3 postProc_check_dSdx.py

echo "=============================="
echo " Finished successfully"
echo "=============================="