vc = 1 # [m/s]
kc = 1 # [m2]
mu_o = 1 # [Pa.s]
L = 1 # [m]

epsilon = 0.01
Pe = 1/epsilon

pc = vc*mu_o*L / (kc*Pe)
print(pc)