#!/bin/bash

# Running it:
# chmod +x build_solvers.sh
# ./build_solvers.sh

# Check if 'solvers' folder exists
if [ ! -d "solvers" ]; then
    echo "Error: 'solvers' folder not found in the current directory."
    exit 1
fi

# Iterate over all directories within 'solvers'
for dir in solvers/*/; do
    if [ -d "$dir" ]; then
        echo "Entering directory: $dir"
        cd "$dir" || { echo "Failed to enter directory: $dir"; exit 1; }
        
        # Run wmake command
        if wmake; then
            echo "wmake successful in $dir"
        else
            echo "wmake failed in $dir"
            exit 1
        fi
        
        # Go back to the starting directory
        cd - > /dev/null || { echo "Failed to return to the base directory"; exit 1; }
    fi
done

echo "All directories processed successfully."
