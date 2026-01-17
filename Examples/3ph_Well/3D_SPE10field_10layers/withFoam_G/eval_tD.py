import numpy as np

def unorm(S_norm, Sr, Swc, Sgr, Sor):
    return S_norm*(1-Swc-Sgr-Sor) + Sr

u = 1.4111e-5

ti = 0
tf = 60000

write_interval = 50

t = np.linspace(ti,tf,int(tf/write_interval + 1))
print(t)

Swc = 0.1
Sgr = 0.0
Sor = 0.1

L = 3.67
phi = 0.25

tD = (u*t) / ((1.0-Swc-Sgr-Sor)*L*phi)

tD_find = 1.0
idx = np.argmin(abs(tD-tD_find))

print(idx, t[idx], tD[idx])


print("INITIAL (I):")
Sw_norm = 0.775
So_norm = 0.225
Sw = unorm(Sw_norm, Swc, Swc, Sgr, Sor)
So = unorm(So_norm, Sor, Swc, Sgr, Sor)
print("Sw = ", Sw, " | So = ", So, " | Sg = ", 1 - Sw - So)
print("\n")
print("INJECTION (J):")
Sw_norm = 0.3125
So_norm = 0.0
Sw = unorm(Sw_norm, Swc, Swc, Sgr, Sor)
So = unorm(So_norm, Sor, Swc, Sgr, Sor)
print("Sw = ", Sw, " | So = ", So, " | Sg = ", 1 - Sw - So)