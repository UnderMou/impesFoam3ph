L = 2000
H = 100
W = 750
phi = 0.07

Vr = L*H*W  # m3

PV = Vr*phi

Qprod_year = 0.05*PV            # PV/year
print(Qprod_year)
year = 365*24*60*60             # s
Qprod_s = Qprod_year / year     # PV/s
print(Qprod_s)

Uinj = 1e-4         # m/s
A = Qprod_s / Uinj  # m2
w = A / H

print(w)

print(Qprod_s/(H*W))
