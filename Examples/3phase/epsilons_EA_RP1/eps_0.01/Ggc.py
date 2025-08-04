import numpy as np 

def kr(S):
    exp = 2
    return S**exp

def pcwo(sw,eps):
    return 5*eps*(2-sw)*(1-sw)

def pcgo(sg,eps):
    return 5*eps*(2-sg)*(1-sg)

eps = 0.02

sw = np.linspace(0.001,0.999,1000)
sg = sw
pcwo_v = pcwo(sw,eps)
pcgo_v = pcgo(sg,eps)

mu_o = 1.0
L = 1.0
vc = 1.0
pc = np.mean(pcwo_v)
kc = 0.5

Pe = vc * mu_o * L / (kc * pc)
eps = 1 / Pe

print('eps\t=\t', eps, '\nkc\t=\t', kc)


# # # Velocities
# # # g = 0
# Sw = 0.613
# Sg = 0.387
# So = 0.0
# krg = kr(Sg)
# krw = kr(Sw)
# kro = kr(So)
# Mg = kc * krg / mu_g
# Mw = kc * krw / mu_w
# Mo = kc * kro / mu_o
# Mt = Mg + Mw + Mo
# fw = Mw/Mt
# fg = Mg/Mt
# # Lt = Mw*rho_w + Mg*rho_g + Mo*rho_o
# # # Lt = 1
# # print("Lt = ", Lt)
# # Lwo = fw*Mo*(rho_w-rho_o)*Kc 
# # Lwg = fw*Mg*(rho_w-rho_g)*Kc
# # Lgo = fg*Mo*(rho_g-rho_o)*Kc 
# # Lgw = fg*Mw*(rho_g-rho_w)*Kc
# # print(Lwo, Lwg, Lgo, Lgw)
# Ug = fg*vc # + (Lwo+Lwg)/(Lt)*(-g)
# Uw = fw*vc # + (Lgo+Lgw)/(Lt)*(-g)
# Uo = vc - Ug - Uw
# print("Ug = ", Ug, "\nUw = ", Uw, "\nUo = ", Uo)
# print("U = ", np.sum(Ug+Uw+Uo))

# ##################################3



# tD  = 0.25

# phi = 0.25
# L   = 1.0
# Swc = 0.0
# Sor = 0.0
# u   = 1.0e-5

# t = tD * (1.0 - Swc - Sor) * L * phi / u 
# print("t : ", t)



# tD  = 0.5

# phi = 0.25
# L   = 1.0
# Swc = 0.0
# Sor = 0.0
# u   = 1.0

# t = tD * (1.0 - Swc - Sor) * L * phi / u 
# print("t : ", t)