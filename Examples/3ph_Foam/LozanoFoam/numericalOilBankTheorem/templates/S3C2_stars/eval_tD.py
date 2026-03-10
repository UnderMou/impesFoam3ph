import numpy as np

def unorm(S_norm, Sr, Swc, Sgr, Sor):
    return S_norm*(1-Swc-Sgr-Sor) + Sr

u = 1.4111e-5

ti = 0
tf = 100000

write_interval = 1000

t = np.linspace(ti,tf,int(tf/write_interval + 1))
print(t)

Swc = 0.0
Sgr = 0.0
Sor = 0.0

L = 1.0
phi = 1.0

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


# GET U inj
def kr(S):
    return S**2

lamb_w = kr(Sw) / 1.0
lamb_o = kr(So) / 4.0
MRF = 1750
lamb_g = kr(1 - Sw - So) / (0.02 * MRF)
lamb_t = lamb_w + lamb_o + lamb_g

fg = lamb_g / lamb_t
fw = lamb_w / lamb_t

uw_inj = u * fw
ug_inj = u * fg
uo_inj = u * (1 - fw - fg)
print("uw_inj = ", uw_inj, " | uo_inj = ", uo_inj ," | ug_inj = ", ug_inj)
print(fw,1-fw-fg,fg)