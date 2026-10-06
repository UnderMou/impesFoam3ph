#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import fluidfoam


# ============================================================
# USER SETTINGS
# ============================================================

case_path = "."           # OpenFOAM case
time_name = "latestTime"  # e.g. "100", "0.5", "latestTime"
field_name = "Sb"

Nx = 51
Ny = 53

# If the 2D plane is x-y:
plane = "xy"

# Show cell borders
show_mesh = False


# ============================================================
# READ OPENFOAM MESH
# ============================================================

print("Reading mesh...")

X, Y, Z = fluidfoam.readmesh(
    case_path,
    structured=True
)

print("Mesh shape:")
print("X =", X.shape)
print("Y =", Y.shape)
print("Z =", Z.shape)


# ============================================================
# READ SATURATION FIELD
# ============================================================

print(f"\nReading {field_name} at {time_name}...")

Sb = fluidfoam.readscalar(
    case_path,
    time_name,
    field_name,
    structured=True
)

print(f"{field_name} shape =", Sb.shape)

print(f"{field_name} min =", np.min(Sb))
print(f"{field_name} max =", np.max(Sb))


# ============================================================
# EXTRACT 2D PLANE
# ============================================================

# For a genuinely 2D OpenFOAM mesh there is normally one cell
# in the third direction.

if plane == "xy":

    X2D = np.squeeze(X)
    Y2D = np.squeeze(Y)
    Sb2D = np.squeeze(Sb)

elif plane == "xz":

    X2D = np.squeeze(X)
    Y2D = np.squeeze(Z)
    Sb2D = np.squeeze(Sb)

elif plane == "yz":

    X2D = np.squeeze(Y)
    Y2D = np.squeeze(Z)
    Sb2D = np.squeeze(Sb)

else:
    raise ValueError("plane must be 'xy', 'xz' or 'yz'")


print("\n2D shape:")
print("X2D  =", X2D.shape)
print("Y2D  =", Y2D.shape)
print("Sb2D =", Sb2D.shape)


# ============================================================
# CHECK Nx x Ny
# ============================================================

if Sb2D.size != Nx * Ny:
    raise ValueError(
        f"Number of cells ({Sb2D.size}) does not match "
        f"Nx*Ny = {Nx}*{Ny} = {Nx*Ny}"
    )


# FluidFoam structured arrays are generally (Nx, Ny)
# We transpose them because matplotlib interprets
# rows as the vertical direction and columns as horizontal.

if Sb2D.shape == (Nx, Ny):

    Xplot = X2D.T
    Yplot = Y2D.T
    Splot = Sb2D.T

elif Sb2D.shape == (Ny, Nx):

    Xplot = X2D
    Yplot = Y2D
    Splot = Sb2D

else:

    # fallback if FluidFoam returned flattened/squeezed data
    Xplot = np.reshape(X2D, (Nx, Ny), order="F").T
    Yplot = np.reshape(Y2D, (Nx, Ny), order="F").T
    Splot = np.reshape(Sb2D, (Nx, Ny), order="F").T


# ============================================================
# REMOVE CELLS FROM LOWER/UPPER Y BOUNDARIES
# ============================================================

notPlot = 1

if notPlot < 0:
    raise ValueError("notPlot must be >= 0")

if 2*notPlot >= Ny:
    raise ValueError(
        f"notPlot={notPlot} removes all cells in y."
    )

if notPlot > 0:

    Xplot = Xplot[notPlot:-notPlot, :]
    Yplot = Yplot[notPlot:-notPlot, :]
    Splot = Splot[notPlot:-notPlot, :]


print("\nPlot dimensions:")
print("Splot shape =", Splot.shape)

print(
    f"Plotting {Nx} x {Ny - 2*notPlot} cells"
)

# ============================================================
# PLOT
# ============================================================

fig, ax = plt.subplots(figsize=(10, 5))

if show_mesh:

    pcm = ax.pcolormesh(
        Xplot,
        Yplot,
        Splot,
        shading="nearest",
        edgecolors="k",
        linewidth=0.15,
        cmap="coolwarm",
        vmin=0.2,
        vmax=0.8
    )

else:

    pcm = ax.pcolormesh(
        Xplot,
        Yplot,
        Splot,
        shading="nearest",
        cmap="coolwarm"#,
        # vmin=0.42,
        # vmax=0.46
    )


# Colorbar
cbar = fig.colorbar(pcm, ax=ax)

cbar.set_label(
    r"$S_b$"
)


# Labels
if plane == "xy":
    ax.set_xlabel("x")
    ax.set_ylabel("y")

elif plane == "xz":
    ax.set_xlabel("x")
    ax.set_ylabel("z")

elif plane == "yz":
    ax.set_xlabel("y")
    ax.set_ylabel("z")


ax.set_title(
    f"{field_name} - time = {time_name}"
)

ax.set_aspect("equal")

plt.tight_layout()
plt.show()