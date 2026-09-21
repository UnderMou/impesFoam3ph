"""Print the OpenFOAM cell indices in the leftmost and rightmost columns."""

import argparse
from pathlib import Path

import fluidfoam
import numpy as np


def column_indices(x, y, z, x_coordinate):
    """Return a column's cell indices, ordered from bottom to top."""
    coordinate_scale = max(float(np.max(np.abs(x))), 1.0)
    tolerance = 16.0 * np.finfo(float).eps * coordinate_scale
    indices = np.flatnonzero(
        np.isclose(x, x_coordinate, rtol=0.0, atol=tolerance)
    )

    # lexsort uses its last key first: y is primary and z breaks ties.
    return indices[np.lexsort((z[indices], y[indices]))]


def get_well_indices(case_directory, ny, nx):
    """Return the cell indices in the leftmost and rightmost columns.

    The mesh is expected to have Ny rows in the y-direction and Nx columns in
    the x-direction.
    """
    x, y, z = fluidfoam.readmesh(
        str(case_directory), structured=False, verbose=False
    )
    x = np.asarray(x).ravel()
    y = np.asarray(y).ravel()
    z = np.asarray(z).ravel()

    if x.size != ny * nx:
        raise ValueError(
            f"Expected {ny * nx} cell centres for Ny={ny}, Nx={nx}, but found {x.size}."
        )

    leftmost = column_indices(x, y, z, float(np.min(x)))
    rightmost = column_indices(x, y, z, float(np.max(x)))

    if leftmost.size != ny or rightmost.size != ny:
        raise ValueError(
            f"Expected {ny} cells in each outer column for Ny={ny}, Nx={nx}, "
            f"but found {leftmost.size} on the left and {rightmost.size} on the right."
        )

    return leftmost, rightmost


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Print the OpenFOAM cell indices in the leftmost and rightmost "
            "columns of a structured mesh."
        )
    )
    parser.add_argument("ny", type=int, help="Number of rows in the y-direction")
    parser.add_argument("nx", type=int, help="Number of columns in the x-direction")
    args = parser.parse_args()

    case_directory = Path(__file__).resolve().parent
    leftmost, rightmost = get_well_indices(case_directory, args.ny, args.nx)

    print("Leftmost column cell indices:")
    print(leftmost)
    print("Rightmost column cell indices:")
    print(rightmost)

    Lx = 2000
    Ly = 100
    Lz = 750
    dx = Lx/args.nx
    dy = Ly/args.ny
    dz = Lz
    dv = dx*dy*dz
    print(1.2208269913749366e-08 * dv * args.ny)

if __name__ == "__main__":
    main()

