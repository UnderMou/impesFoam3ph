#!/bin/bash

echo "Starting compilation..."

set -e

dirs=(
    "./boundaryConditions/2ph"
    "./boundaryConditions/3ph"
    "./capillaryPressureModels/2ph"
    "./capillaryPressureModels/3ph"
    "./foamModels/2ph"
    "./foamModels/3ph"
    "./isothermModels/2ph"
    "./isothermModels/3ph"
    "./relativePermeabilityModels/2ph"
    "./relativePermeabilityModels/3ph"
    "./surfactantTransportModels/2ph"
    "./surfactantTransportModels/3ph"
    "./wellModels/2ph"
    "./wellModels/3ph"
)

for module_dir in "${dirs[@]}"; do
    if [ -d "$module_dir" ] && [ -d "$module_dir/Make" ]; then
        echo "--------------------------------"
        echo ""

        cd "$module_dir" || exit

        pwd
        echo ""

        wclean
        wmake libso .

        cd - > /dev/null
    else
        echo "Skipping $module_dir (directory or Make folder not found)"
    fi
done
echo "--------------------------------"

for dir in solvers/*; do
    if [ -d "$dir" ]; then
        echo ""
        
        cd "$dir" || exit

        pwd
        echo ""

        wclean
        wmake

        echo "--------------------------------"

        cd - > /dev/null
    fi
done

echo "Compilation finished."