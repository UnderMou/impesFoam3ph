import numpy as np
import csv
import ternary
import matplotlib.pyplot as plt

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
    Sw_min=0.05, Sw_max=1.0,
    So_min=0.002, So_max=1.0,
    Sg_min=0.0, Sg_max=1.0
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

def injection_points(k,
    Sw_min=0.0, Sw_max=0.7,
    So_min=0.0, So_max=0.0,
    Sg_min=0.3, Sg_max=1.0):

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
    

    # Graph
    scale = 1
    figure, tax = ternary.figure(scale=scale)
    tax.gridlines(multiple=0.1, linewidth=0.5, color="blue")
    tax.boundary(linewidth=2.0)
    fontsize = 14
    tax.right_corner_label("W", fontsize=fontsize, color="blue", fontweight="bold")
    tax.top_corner_label("O", fontsize=fontsize, color="black", fontweight="bold")
    tax.left_corner_label("G", fontsize=fontsize, color="green", fontweight="bold")

    # INITIAL CONDITIONS
    k = 25
    # Sw, So, Sg = ternary_grid(k)
    Sw, So, Sg = ternary_grid_bounded(k)
    points = [(float(w)*scale, float(o)*scale, float(g)*scale) for w,o,g in zip(Sw,So,Sg)]
    tax.scatter(points, marker='s', color='dimgray', label="Initial")
    print("nR points: ", len(points))

    with open("initial_conditions.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sw", "so", "sg"])  # header
        writer.writerows(points)


    # INJECTION CONDITIONS
    k = 10
    points = injection_points(k)
    tax.scatter(points, marker='s', color='red', label="Injection")
    
    tax.ticks(axis='lbr', linewidth=1, offset=0.025, multiple=0.1, tick_formats="%.1f")
    tax.get_axes().axis('off')
    tax.clear_matplotlib_ticks()
    print("nL points: ", len(points))
    print(points)

    with open("injection_conditions.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sw", "so", "sg"])  # header
        writer.writerows(points)

    
    plt.legend()
    tax.show()