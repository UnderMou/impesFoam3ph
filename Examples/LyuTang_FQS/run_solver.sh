#!/bin/bash

source /opt/openfoam9/etc/bashrc

if [ -z "$1" ]; then
    echo "Error: missing base case directory"
    exit 1
fi

if [ ! -d "$1" ]; then
    echo "Error: case directory does not exist"
    exit 1
fi

# # Check if the destination directory is provided
# if [ -z "$2" ]; then
# echo "Error: missing new case directory"
# exit 1
# fi

# # Create the new directory and copy only the contents of the base case into it
# mkdir -p "$2"
# cp -a "$1/." "$2"

# Move to the new case directory and run the solver
cd "$1"
impesFoam3ph_v1