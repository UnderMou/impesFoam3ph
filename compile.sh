#!/bin/bash

echo "Starting compilation..."

# Find all directories containing a Make folder,
# but skip solvers and Examples
# find . \
#     -path "./solvers" -prune -o \
#     -path "./Examples" -prune -o \
#     -type d -name Make -print | while read dir; do

#     module_dir=$(dirname "$dir")

#     echo "--------------------------------"
#     echo ""
#     cd "$module_dir" || exit

#     # Print current directory
#     pwd
#     echo ""

#     # Compile
#     wclean
#     wmake libso .

#     cd - > /dev/null
# done

set -e

dirs=(
    "./boundaryConditions/2ph"
    "./boundaryConditions/3ph"
    "./capillaryPressureModels/2ph"
    "./foamModels/2ph"
    "./foamModels/3ph"
    "./isothermModels"
    "./relativePermeabilityModels/2ph"
    "./relativePermeabilityModels/3ph"
    "./surfactantTransportModels"
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