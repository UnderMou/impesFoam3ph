import numpy as np
import matplotlib.pyplot as plt

def pcgo(sw,sg):
    return 5*(2-sg)*(1-sg)

def pcwo(sw,sg):
    return 5*(2-sw)*(1-sw)

def kri(Si):
    return Si**2

def DpwoDSw(sw,sg):
    return 5*(2*sw-3)

def DpwoDSg(sw,sg):
    return 0

def DpgoDSw(sw,sg):
    return 0

def DpgoDSg(sw,sg):
    return 5*(2*sg-3)

mu_g = 0.3
mu_w = 0.5
mu_o = 1.0

Sw = np.linspace(0.001,0.999,100)
Sg = Sw

fig, ax = plt.subplots(1,2,figsize=(7,5))
ax[0].plot(Sw, pcwo(Sw,Sg))
ax[1].plot(Sg, pcgo(Sw,Sg))
plt.show()

Sw_grid, Sg_grid = np.meshgrid(Sw, Sg)
S = np.column_stack((Sw_grid.ravel(), Sg_grid.ravel()))
for i in range(S.shape[0]):
    sw = S[i][0]
    sg = S[i][1]
    lambda_w = kri(sw)
    lambda_g = kri(sg)

    Q = np.array([[lambda_w, 0       ],
                [0       , lambda_g]])
    P = np.array([[DpwoDSw(sw,sg), DpwoDSg(sw,sg)],
                [DpgoDSw(sw,sg), DpgoDSg(sw,sg)]])
    B = Q @ P

    det_B = np.linalg.det(B)

    if det_B <= 0:
        print(sw, sg)