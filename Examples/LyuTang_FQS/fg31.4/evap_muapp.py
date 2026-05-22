import numpy as np

L = 0.4 # [m]
Pin = 6.89e6 # [Pa]
Pout = 5e6 # [Pa]
K = 2.783123706e-12 # [m2]
ut = 1.61572e-5 # [m/s]

gradP = (Pin - Pout) / L # [Pa/m]
muapp = (K * gradP) / (ut) # [Pa.s]

# muapp *= 1e3 # [cP]
print(muapp)