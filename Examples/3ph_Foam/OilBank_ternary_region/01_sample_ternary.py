import numpy as np
import csv
import ternary
import matplotlib.pyplot as plt
import json

def ternary_grid(k):
    Sw, So, Sg = [], [], []

    for i in range(k+1):
        for j in range(k+1-i):
            l = k - i - j

            Sw.append(i/k)
            So.append(j/k)
            Sg.append(l/k)

    return np.array(Sw), np.array(So), np.array(Sg)

def ternary_grid_bounded(
    k,
    Sw_min,
    Sw_max,
    So_min,
    So_max,
    Sg_min,
    Sg_max
):
    Sw, So, Sg = [], [], []

    for i in range(k+1):
        for j in range(k+1-i):
            l = k - i - j

            sw = i / k
            so = j / k
            sg = l / k

            if (
                Sw_min <= sw <= Sw_max and
                So_min <= so <= So_max and
                Sg_min <= sg <= Sg_max
            ):
                Sw.append(sw)
                So.append(so)
                Sg.append(sg)

    return np.array(Sw), np.array(So), np.array(Sg)

def injection_points(
    k,
    Sw_min,
    Sw_max,
    So_min,
    So_max,
    Sg_min,
    Sg_max
):

    pts = []

    for i in range(k+1):
        for j in range(k+1-i):
            l = k - i - j

            Sw = i/k
            So = j/k
            Sg = l/k

            if (
                Sw_min <= Sw <= Sw_max and
                So_min <= So <= So_max and
                Sg_min <= Sg <= Sg_max
            ):
                pts.append((Sw, So, Sg))

    return pts


if __name__ == "__main__":

    # User defined ===========
    case_name = "CTscan_moreInj"
    # ========================

    path2save = f"casesPrep/{case_name}/"
    
    with open("initials_range.json", "r") as f:
        initials_range = json.load(f)
    
    with open("injections_range.json", "r") as f:
        injections_range = json.load(f)

    initials_range = initials_range[case_name]
    injections_range = injections_range[case_name]

    # Graph
    scale = 1
    figure, tax = ternary.figure(scale=scale)
    tax.gridlines(multiple=0.1, linewidth=0.5, color="blue")
    tax.boundary(linewidth=2.0)
    
    # INITIAL CONDITIONS
    Sw, So, Sg = ternary_grid_bounded(**initials_range)
    points = [(float(w)*scale, float(o)*scale, float(g)*scale) for w,o,g in zip(Sw,So,Sg)]
    tax.scatter(points, marker='s', color='dimgray', label="Initial")
    nR = len(points)
    print("nR points: ", nR)

    with open(f"{path2save}/initial_conditions_{case_name}.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sw", "so", "sg"])  # header
        writer.writerows(points)


    # INJECTION CONDITIONS
    points = injection_points(**injections_range)
    tax.scatter(points, marker='s', color='red', label="Injection")
    
    tax.ticks(axis='lbr', linewidth=1, offset=0.025, multiple=0.1, tick_formats="%.1f")
    tax.get_axes().axis('off')
    tax.clear_matplotlib_ticks()
    nL = len(points)
    print("nL points: ", nL)
    print(points)

    print("total simulations: ", nR*nL)

    with open(f"{path2save}/injection_conditions_{case_name}.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sw", "so", "sg"])  # header
        writer.writerows(points)


    # PLOT
    fontsize = 14
    tax.right_corner_label("W", fontsize=fontsize)
    tax.top_corner_label("O", fontsize=fontsize)
    tax.left_corner_label("G", fontsize=fontsize)
    
    tax.legend()
    tax.savefig(f"{path2save}/samples_{case_name}.png", dpi=300)
    tax.show()
