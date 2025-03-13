import numpy as np

def unorm(S_norm, Sr, Swc, Sgr, Sor):
    return S_norm*(1-Swc-Sgr-Sor) + Sr

# L = 0.4
# D = 0.04
# A = np.pi * np.power(D/2,2)
# V = A * L

u = 1.61572e-5

ti = 0
tf = 10000

write_interval = 25

t = np.linspace(ti,tf,int(tf/write_interval + 1))
print(t)

Swc = 0.197
Sgr = 0.013
Sor = 0.103

L = 0.4
phi = 0.22

tD = (u*t) / ((1.0-Swc-Sgr-Sor)*L*phi)
# tD = (u*t*A) / ((1.0-Swc-Sgr-Sor)*V*phi)

print(tD)

tD_find = 0.36
idx = np.argmin(abs(tD-tD_find))

print(idx, t[idx], tD[idx])


print("INITIAL (I):")
Sw_norm = 0.54
So_norm = 0.46
Sw = unorm(Sw_norm, Swc, Swc, Sgr, Sor)
So = unorm(So_norm, Sor, Swc, Sgr, Sor)
print("Sw = ", Sw, " | So = ", So, " | Sg = ", 1 - Sw - So)
print("\n")
# print("INJECTION (J):")
# Sw_norm = 0.2
# So_norm = 0.0
# Sw = unorm(Sw_norm, Swc, Swc, Sgr, Sor)
# So = unorm(So_norm, Sor, Swc, Sgr, Sor)
# print("Sw = ", Sw, " | So = ", So, " | Sg = ", 1 - Sw - So)