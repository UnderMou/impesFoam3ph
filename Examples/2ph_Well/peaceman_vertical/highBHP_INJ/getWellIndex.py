"""Print the OpenFOAM cell indices in the leftmost and rightmost columns."""

from pathlib import Path

import fluidfoam
import numpy as np


NUMBER_OF_ROWS = 51


def column_indices(x, y, z, x_coordinate):
    """Return a column's cell indices, ordered from bottom to top."""
    coordinate_scale = max(float(np.max(np.abs(x))), 1.0)
    tolerance = 16.0 * np.finfo(float).eps * coordinate_scale
    indices = np.flatnonzero(
        np.isclose(x, x_coordinate, rtol=0.0, atol=tolerance)
    )

    # lexsort uses its last key first: y is primary and z breaks ties.
    return indices[np.lexsort((z[indices], y[indices]))]


def main():
    case_directory = Path(__file__).resolve().parent

    # fluidfoam resolves this case path to constant/polyMesh and returns the
    # cell centres in OpenFOAM cell-index order when structured=False.
    x, y, z = fluidfoam.readmesh(
        str(case_directory), structured=False, verbose=False
    )
    x = np.asarray(x).ravel()
    y = np.asarray(y).ravel()
    z = np.asarray(z).ravel()

    leftmost = column_indices(x, y, z, float(np.min(x)))
    rightmost = column_indices(x, y, z, float(np.max(x)))

    if leftmost.size != NUMBER_OF_ROWS or rightmost.size != NUMBER_OF_ROWS:
        raise ValueError(
            "Expected 51 cells in each outer column, but found "
            f"{leftmost.size} on the left and {rightmost.size} on the right."
        )

    print("Leftmost column cell indices:")
    print(leftmost)
    print("Rightmost column cell indices:")
    print(rightmost)


if __name__ == "__main__":
    main()