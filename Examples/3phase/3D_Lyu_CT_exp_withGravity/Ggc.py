import numpy as np 

def kr(S):
    exp = 2
    return S**exp


# Ggc     = 0.0055
g       = 9.81
rho_o   = 0.7
rho_w   = 1.0
rho_g   = 5.76e-2
rho_wgc = rho_w - rho_g
vc      = 1.1310039999999998e-05 + 4.84716e-06
mu_o    = 9.0e-2
Kc      = 2.783123706e-12

Ggc = Kc * rho_wgc * g / ( vc * mu_o)
print("Ggc = ", Ggc)
exit()


# Kc = mu_o*Ggc*vc / (rho_wgc*g)
mu_w = 0.7e-3
mu_g = 1.5e-5
# print(f"mu_o : {mu_o:.2e}\nmu_w : {mu_w:.2e}")


# vc = Kc * rho_wgc * g / (Ggc * mu_o)
# print(f"vc = {vc:.2e}")



# PECLET

eps = 0.001
Pe = 1.0 / eps
# TODO : tem q implementar pressoes capilares


# Velocities
# g = 0
Sw = 0.721
Sg = 0.279
So = 0.0
krg = kr(Sg)
krw = kr(Sw)
kro = kr(So)
Mg = Kc * krg / mu_g
Mw = Kc * krw / mu_w
Mo = Kc * kro / mu_o
Mt = Mg + Mw + Mo
fw = Mw/Mt
fg = Mg/Mt

# Lt = Mw*rho_w + Mg*rho_g + Mo*rho_o
# Lt = 1
# print("Lt = ", Lt)

# Lwo = fw*Mo*(rho_w-rho_o)*Kc 
# Lwg = fw*Mg*(rho_w-rho_g)*Kc
# Lgo = fg*Mo*(rho_g-rho_o)*Kc 
# Lgw = fg*Mw*(rho_g-rho_w)*Kc
# print(Lwo, Lwg, Lgo, Lgw)
# Ug = fg*vc + (Lwo+Lwg)*(-g)
# Uw = fw*vc + (Lgo+Lgw)*(-g)
# Uo = vc - Ug - Uw

Lw = Mw*rho_w
Lg = Mg*rho_g
Uw = fw*vc 
Ug = fg*vc 
Uo = vc - Ug - Uw

print("Ug = ", Ug, "\nUw = ", Uw, "\nUo = ", Uo)
print("U = ", np.sum(Ug+Uw+Uo))



##################################3



# tD  = 0.25

# phi = 0.25
# L   = 1.0
# Swc = 0.0
# Sor = 0.0
# u   = 1.0e-5

# t = tD * (1.0 - Swc - Sor) * L * phi / u 
# print("t : ", t)



tD  = 0.5

phi = 0.25
L   = 1.0
Swc = 0.0
Sor = 0.0
u   = 1.0

t = tD * (1.0 - Swc - Sor) * L * phi / u 
print("t : ", t)

 
print("Ggc = ", Kc * rho_wgc * g / (vc * mu_o))