import meshio
import numpy as np
from itertools import combinations

mesh = meshio.read("UNISIM-I-D_optionA.vtu")

hexes = mesh.cells_dict["hexahedron"]
points = mesh.points

bad_cells = [7174, 7911, 7422, 7505, 7148]

for cell_id in bad_cells:

    conn = hexes[cell_id]
    xyz = points[conn]

    print("\n" + "="*70)
    print(f"CELL {cell_id}")
    print("="*70)

    print("\nConnectivity:")
    print(conn)

    print("\nCoordinates:")
    for local_id, (point_id, coord) in enumerate(zip(conn, xyz)):
        print(
            f"local {local_id}: "
            f"point {point_id:6d} -> "
            f"{coord[0]:14.6f} "
            f"{coord[1]:14.6f} "
            f"{coord[2]:14.6f}"
        )

    # Check every pair of vertices
    print("\nVery close / coincident vertices:")

    found = False

    for i, j in combinations(range(8), 2):

        d = np.linalg.norm(xyz[i] - xyz[j])

        if d < 1e-6:
            print(
                f"  local {i} - local {j}: "
                f"distance = {d:.12e}"
            )
            found = True

    if not found:
        print("  None")