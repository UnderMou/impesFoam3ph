import numpy as np 

def kr(S):
    exp = 2
    return S**exp

# Ggc     = 2
# g       = 9.81
# rho_w   = 1.0
# rho_o   = 0.7
# rho_woc = rho_w - rho_o
# vc      = 1.0e-5
# # mu_o    = 1.0e-3
# Kc      = 1.0e-13

# # Kc = Ggc * vc * mu_o / (rho_woc * g)
# # print(f"{Kc:.2e}")

# mu_o = Kc*rho_woc*g/(Ggc*vc)
# mu_w = mu_o / 5
# print(f"mu_o : {mu_o:.2e}\nmu_w : {mu_w:.2e}")

# Ggc     = 0.0055
Ggc     = 0.085
g       = 9.81
rho_o   = 0.7
rho_w   = 1.0
rho_g   = 5.76e-2
rho_wgc = rho_w - rho_g
vc      = 1.0
mu_o    = 1.0
# Kc      = 1.0

Kc = Ggc * vc * mu_o / (rho_wgc * g)
# Kc = 1e-1
print(f"Kc = {Kc:.2e}")



# Kc = mu_o*Ggc*vc / (rho_wgc*g)
mu_w = 0.5
mu_g = 0.3
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