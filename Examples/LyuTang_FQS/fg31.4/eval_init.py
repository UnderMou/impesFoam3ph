import numpy as np

def unorm(S_norm, Sr, Swc, Sgr, Sor):
    return S_norm*(1-Swc-Sgr-Sor) + Sr

Swc = 0.135 #0.197
Sgr = 0.2   #0.013
Sor = 0.1   #0.103

Swi_norm = 0.54
Soi_norm = 0.46
Sgi_norm = 1.0 - Swi_norm - Soi_norm
Swi = unorm(Swi_norm, Swc, Swc, Sgr, Sor)
Soi = unorm(Soi_norm, Sor, Swc, Sgr, Sor)
Sgi = unorm(Sgi_norm, Sgr, Swc, Sgr, Sor)

print("Swi: ", Swi, "\nSoi: ", Soi, "\nSgi: ", Sgi)