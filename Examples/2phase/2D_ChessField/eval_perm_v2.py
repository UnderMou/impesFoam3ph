import pandas as pd
import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import cv2

# ts = np.linspace(0,150000,601, endpoint=True, dtype=int)
# print(ts.tolist())
nx = 200
ny = 60

# Mapping the permeability field
min_val, max_val = 1e-13, 4e-13
image_path = "perm_quad.png"  # Replace with the path to your image
image = cv2.imread(image_path)
gray_image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
gray_matrix = np.array(gray_image)
inverted_matrix = gray_matrix[::-1, :]
gray_matrix = inverted_matrix
print(gray_matrix.shape)
normalized_gray = gray_matrix / 255.0  # Values from 0 (black) to 1 (white)
real_value_matrix = normalized_gray * (max_val - min_val) + min_val

new_shape = (nx, ny)  # Target shape: (width, height)
perm = cv2.resize(real_value_matrix, new_shape, interpolation=cv2.INTER_LINEAR)

# WRITE OPENFOAM-LIKE FILE

# Path to write the OpenFOAM-like file
file_path = 'K'

# Open file and write in the OpenFOAM format
with open(file_path, 'w') as f:
    f.write(
        '''/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  9
     \\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       volScalarField;
    location    "0";
    object      K;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar> 
'''
    )
    f.write(f"{nx*ny}\n(\n")
    for row in perm:
        for value in row:
            f.write(f"{value}\n")  # Write each perm value in scientific notation
    f.write(")\n")

    f.write(
        ''';
boundaryField
{
    inlet
    {
        type            zeroGradient;
    }
    outlet
    {
        type            zeroGradient;
    }
    walls
    {
        type            zeroGradient;
    }
    emptyWalls
    {
        type            empty;
    }
}'''
    )

print(f"File written to {file_path}")