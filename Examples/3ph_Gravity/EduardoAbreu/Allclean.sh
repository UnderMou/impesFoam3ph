#!/bin/bash

echo "=============================="
echo " Cleaning case directories"
echo "=============================="

# Remove files (ignore if they don't exist)
rm -f *.pdf

# List of case directories
dirs=(CaseA CaseB CaseB_noG)

# Store root directory
ROOT_DIR=$(pwd)

for dir in "${dirs[@]}"; do
    cd "$ROOT_DIR/$dir" || exit 1

    echo "cleaning directory: $dir"
    foamListTimes -rm 
    echo ""

done

echo "=============================="
echo " Cleanup finished"
echo "=============================="