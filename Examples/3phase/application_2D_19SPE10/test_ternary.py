import numpy as np
import matplotlib.pyplot as plt
import ternary
import matplotlib.tri as tri


# --- Triangle corners ---
corner_So = np.array([0.0, 0.0])             
corner_Sg = np.array([1.0, 0.0])             
corner_Sw = np.array([0.5, np.sqrt(3)/2])    

# Resolution
N = 10
So_vals = np.linspace(0, 1, N)
Sg_vals = np.linspace(0, 1, N)

points = []
So_list, Sg_list, Sw_list = [], [], []

for So in So_vals:
    for Sg in Sg_vals:
        Sw = 1 - So - Sg
        if Sw >= 0:
            x = So * corner_So[0] + Sg * corner_Sg[0] + Sw * corner_Sw[0]
            y = So * corner_So[1] + Sg * corner_Sg[1] + Sw * corner_Sw[1]
            points.append([x, y])
            So_list.append(So)
            Sg_list.append(Sg)
            Sw_list.append(Sw)

points = np.array(points)
So_arr, Sg_arr, Sw_arr = np.array(So_list), np.array(Sg_list), np.array(Sw_list)

# ---- Clip + rescale ----
# # OPTION 1
# So_min, So_max = 0., 0.44
# Sg_min, Sg_max = 0., 0.61
# Sw_min, Sw_max = 0., 0.72
# OPTION 2
So_min, So_max = 0., 0.28
Sg_min, Sg_max = 0., 0.18
Sw_min, Sw_max = 0., 0.72
# # OPTION 3
# So_min, So_max = 0.1, 0.28
# Sg_min, Sg_max = 0., 0.18
# Sw_min, Sw_max = 0.57, 0.74

def clip_and_rescale(arr, vmin, vmax):
    arr_clipped = np.clip(arr, vmin, vmax)
    return (arr_clipped - vmin) / (vmax - vmin + 1e-12)

So_scaled = clip_and_rescale(So_arr, So_min, So_max)
Sg_scaled = clip_and_rescale(Sg_arr, Sg_min, Sg_max)
Sw_scaled = clip_and_rescale(Sw_arr, Sw_min, Sw_max)

triang = tri.Triangulation(points[:,0], points[:,1])

# --- Plot triangle ---
fig, ax = plt.subplots(figsize=(4,4))
ax.set_aspect("equal")

vertex_colors = np.vstack([So_scaled, Sg_scaled, Sw_scaled]).T
for t in triang.triangles:
    poly = points[t]
    c = vertex_colors[t].mean(axis=0)
    ax.fill(poly[:,0], poly[:,1], color=c, edgecolor='none')

# Triangle border
triangle = np.array([corner_So, corner_Sg, corner_Sw, corner_So])
ax.plot(triangle[:,0], triangle[:,1], 'k-', lw=1.5)

# Labels
ax.text(*corner_So, "$S_o$", color="red", ha="right", va="top", fontsize=22)
ax.text(*corner_Sg, "$S_g$", color="green", ha="left", va="top", fontsize=22)
ax.text(*corner_Sw, "$S_w$", color="blue", ha="center", va="bottom", fontsize=22)

# --- Add ticks on each side ---
num_ticks = 5
offset = 0.03
for i in range(1, num_ticks):
    frac = i / num_ticks
    
    # Side So-Sg
    x_tick = corner_So[0] + frac * (corner_Sg[0] - corner_So[0])
    y_tick = corner_So[1] + frac * (corner_Sg[1] - corner_So[1])
    So_val = So_min + (1-frac)*(So_max - So_min)
    Sg_val = Sg_min + frac*(Sg_max - Sg_min)
    ax.text(x_tick, y_tick - 2.0*offset, f"{So_val:.2f}", color="red", fontsize=16, ha="center", va="top", rotation=45, rotation_mode="anchor")
    
    # Side Sg-Sw
    x_tick = corner_Sg[0] + frac * (corner_Sw[0] - corner_Sg[0])
    y_tick = corner_Sg[1] + frac * (corner_Sw[1] - corner_Sg[1])
    Sg_val = Sg_min + (1-frac)*(Sg_max - Sg_min)
    Sw_val = Sw_min + frac*(Sw_max - Sw_min)
    ax.text(x_tick + offset, y_tick, f"{Sg_val:.2f}", color="green", fontsize=16, ha="left", va="center")
    
    # Side Sw-So
    x_tick = corner_Sw[0] + frac * (corner_So[0] - corner_Sw[0])
    y_tick = corner_Sw[1] + frac * (corner_So[1] - corner_Sw[1])
    Sw_val = Sw_min + (1-frac)*(Sw_max - Sw_min)
    So_val = So_min + frac*(So_max - So_min)
    ax.text(x_tick - offset, y_tick, f"{Sw_val:.2f}", color="blue", fontsize=16, ha="right", va="center")

ax.set_axis_off()
plt.savefig('ternary.png',dpi=500)
plt.show()
