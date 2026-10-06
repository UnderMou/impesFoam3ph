#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from fluidfoam import readmesh, readscalar
import scienceplots
import matplotlib

plt.style.use('science')
matplotlib.rcParams.update({'font.size': 22})


# ============================================================
# USER SETTINGS
# ============================================================

case_dir = "."
time_name = "2.4e+07" # "2.4e+07" # "9e+06"
field_name = "Sb"

# Coordinates of the debug cell (cell 5)
cellDebug = np.array([15.3401, 0.217027, 5.0])

# Vertical interval
y_min = 0.0
y_max = 143.0

# Tolerances used to identify cells in the same x-z column
tol_x = 1e-4
tol_z = 1e-4

# Optional wall patch
read_wall = True
wall_patch = "walls"


# ============================================================
# READ INTERNAL MESH
# ============================================================

x, y, z = readmesh(
    case_dir,
    structured=False,
    verbose=False
)

x = np.asarray(x).ravel()
y = np.asarray(y).ravel()
z = np.asarray(z).ravel()


# ============================================================
# READ INTERNAL Sb
# ============================================================

Sb = readscalar(
    case_dir,
    time_name,
    field_name,
    structured=False,
    verbose=False
)

Sb = np.asarray(Sb).ravel()

if not (len(x) == len(y) == len(z) == len(Sb)):
    raise RuntimeError(
        f"Inconsistent sizes: x={len(x)}, y={len(y)}, "
        f"z={len(z)}, Sb={len(Sb)}"
    )

print(f"Number of cells: {len(Sb)}")


# ============================================================
# FIND CELL CLOSEST TO cellDebug
# ============================================================

distance = np.sqrt(
    (x - cellDebug[0])**2
    + (y - cellDebug[1])**2
    + (z - cellDebug[2])**2
)

debugCell = np.argmin(distance)

print("\n==============================")
print("DEBUG CELL")
print("==============================")
print(f"index    = {debugCell}")
print(f"C        = ({x[debugCell]:.8f}, "
      f"{y[debugCell]:.8f}, {z[debugCell]:.8f})")
print(f"Sb       = {Sb[debugCell]:.10f}")
print(f"distance = {distance[debugCell]:.6e}")


# ============================================================
# SELECT VERTICAL COLUMN AT SAME x,z
# ============================================================

x0 = x[debugCell]
z0 = z[debugCell]

mask = (
    np.isclose(x, x0, atol=tol_x, rtol=0.0)
    & np.isclose(z, z0, atol=tol_z, rtol=0.0)
    & (y >= y_min)
    & (y <= y_max)
)

column_indices = np.where(mask)[0]

if len(column_indices) == 0:
    raise RuntimeError(
        "No cells were found in the vertical column. "
        "Increase tol_x and/or tol_z."
    )

# Sort by y
column_indices = column_indices[np.argsort(y[column_indices])]

y_line = y[column_indices]
Sb_line = Sb[column_indices]


# ============================================================
# PRINT COLUMN VALUES
# ============================================================

print("\n==============================")
print("VERTICAL COLUMN")
print("==============================")
print(
    f"{'cell':>8} "
    f"{'x':>14} "
    f"{'y':>14} "
    f"{'z':>14} "
    f"{'Sb':>14}"
)

for celli in column_indices:
    print(
        f"{celli:8d} "
        f"{x[celli]:14.6f} "
        f"{y[celli]:14.6f} "
        f"{z[celli]:14.6f} "
        f"{Sb[celli]:14.8f}"
    )


# ============================================================
# OPTIONAL: READ WALL PATCH
# ============================================================

wall_y = None
wall_Sb = None

if read_wall:
    try:
        xw, yw, zw = readmesh(
            case_dir,
            structured=False,
            boundary=wall_patch,
            verbose=False
        )

        Sbw = readscalar(
            case_dir,
            time_name,
            field_name,
            structured=False,
            boundary=wall_patch,
            verbose=False
        )

        xw = np.asarray(xw).ravel()
        yw = np.asarray(yw).ravel()
        zw = np.asarray(zw).ravel()
        Sbw = np.asarray(Sbw).ravel()

        wallPoint = np.array([x0, 0.0, z0])

        wall_distance = np.sqrt(
            (xw - wallPoint[0])**2
            + (yw - wallPoint[1])**2
            + (zw - wallPoint[2])**2
        )

        wallFace = np.argmin(wall_distance)

        wall_y = yw[wallFace]
        wall_Sb = Sbw[wallFace]

        print("\n==============================")
        print("WALL FACE")
        print("==============================")
        print(f"patch    = {wall_patch}")
        print(f"face     = {wallFace}")
        print(
            f"Cf       = ({xw[wallFace]:.8f}, "
            f"{yw[wallFace]:.8f}, {zw[wallFace]:.8f})"
        )
        print(f"Sb       = {Sbw[wallFace]:.10f}")
        print(f"distance = {wall_distance[wallFace]:.6e}")

    except Exception as exc:
        print("\nCould not read wall patch:")
        print(exc)


# ============================================================
# PLOT
# ============================================================

plt.figure(figsize=(8, 7))

plt.plot(
    y_line[1:-1],
    Sb_line[1:-1],
    # y_line,
    # Sb_line,
    "-o",
    linewidth=2,
    markersize=4,
    label="OpenFOAM cell-centre values"
)

plt.scatter(
    [y_line[0], y_line[-1]],
    [Sb_line[0], Sb_line[-1]],
    s=40,
    c="r",
    label="Removed"
)

# # Highlight debug cell
# plt.scatter(
#     y[debugCell],
#     Sb[debugCell],
#     s=100,
#     marker="x",
#     label=f"Debug cell {debugCell}",
#     zorder=10
# )

# # Highlight wall value if successfully read
# if wall_y is not None:
#     plt.scatter(
#         wall_y,
#         wall_Sb,
#         s=100,
#         marker="s",
#         label=f"Patch '{wall_patch}'",
#         zorder=10
#     )

plt.xlabel(r"$y$")
plt.ylabel(r"$S_b$")
plt.xlim(y_min, y_max)
plt.grid(True)
plt.ylim([0.4,0.45])
plt.legend()
plt.tight_layout()

plt.savefig("Sb_vertical_profile.png", dpi=300)
plt.show()
