import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from fluidfoam import readmesh
from fluidfoam import readvector, readscalar
from fluidfoam.readof import getVolumes

import os

sol = '.'
x, y, z = readmesh(sol, structured=True)
nx, ny, nz = x.shape
print("Nx = ", nx, "Ny = ", ny, "Nz = ", nz)


# Mesh
mybox = ((0, 0, 0),(360, 100, 10))
centroidList_box, vol_box = getVolumes( path =  sol, box = mybox)
vol_box = vol_box.reshape((nx, ny, nz))

times = sorted([d for d in os.listdir(sol) if d.isdigit()], key=int)
times = times[1:]

# Figure setup
fig, ax = plt.subplots(figsize=(10,4))

# levels = np.arange(-5e-8, 5e-8, 1e-8)
# levels = np.arange(0, 1.01, 0.001)
levels = np.linspace(-1e-8,1.01,100)
print(levels)

# Initial field (first timestep)
qt = readscalar(sol, times[0], 'qt', structured=True)
prod_wells_idx = qt < 0

Fa0 = readscalar(sol, str(times[0]), 'Fa', structured=True)
Fb0 = readscalar(sol, str(times[0]), 'Fb', structured=True)
Fc0 = 1 - Fa0 - Fb0

contour = ax.contourf(
    x[:, :, nz//2],
    y[:, :, nz//2],
    Fc0[:, :, nz//2],
    levels=levels,
    cmap='coolwarm'
)

cbar = plt.colorbar(contour, ax=ax)
cbar.set_label(r'$f_o$')

ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')

qo = []

# Update function
def update(i):
    ax.clear()
    
    timename = times[i]
    print(f"Processing time: {timename}")
    
    Fa = readscalar(sol, str(times[i]), 'Fa', structured=True)
    Fb = readscalar(sol, str(times[i]), 'Fb', structured=True)
    Fc = 1 - Fa - Fb

    qo.append(-np.sum(Fc[prod_wells_idx] * qt[prod_wells_idx] * vol_box[prod_wells_idx]))
    # print(qo[-1])
    # exit()
    
    contour = ax.contourf(
        x[:, :, nz//2],
        y[:, :, nz//2],
        Fc[:, :, nz//2],
        levels=levels,
        cmap='coolwarm'
    )
    
    ax.set_title(f'Time = {timename}')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    
    return []

# Create animation
anim = FuncAnimation(
    fig,
    update,
    frames=len(times),
    interval=200,   # milliseconds between frames
    blit=False
)

# Save animation
anim.save('animation.mp4', dpi=150, fps=5)
plt.close()

times = np.array(times, dtype=int)
dt = times[1:] - times[:-1]
dt = np.concatenate(([0],dt))
Qo_t = qo[1:] * dt
Qo_cuml = [np.sum(Qo_t[:i]) for i in range(len(Qo_t))]

print(times)
print(Qo_cuml)


plt.plot(times, Qo_cuml)
plt.xlabel(r"$t$ [s]")
plt.ylabel(r"Cumulative oil production [m$^3$]")
plt.grid()
plt.show()

# Or show interactively
# plt.show()

# #tuple of box's dimension: ((xmin, ymin, zmin), (xmax, ymax, zmax))
# mybox = ((0, 0, 0),(360, 100, 10))

# #getVolumes function returns arrays containing the centroids and volume of the
# #cells inside boxes A and B
# centroidList_box,vol_box = getVolumes( path =  sol, box = mybox)

# print(np.mean(vol_box),np.min(vol_box),np.max(vol_box))

# vol_box = vol_box.reshape((nx, ny, nz))
# lx, ly, lz = 360, 100, 10
# dx = lx / nx
# dy = ly / ny
# dz = lz / nz
# dv = dx * dy * dz
# print(dv)


# im = plt.contourf(
#     x[:, :, nz//2],
#     y[:, :, nz//2],
#     vol_box[:, :, nz//2],
#     cmap='coolwarm'
# )

# plt.colorbar(im, label=r'$F_a$')

# # Setting axis labels
# plt.xlabel('x (m)')
# plt.ylabel('y (m)')

# plt.show()