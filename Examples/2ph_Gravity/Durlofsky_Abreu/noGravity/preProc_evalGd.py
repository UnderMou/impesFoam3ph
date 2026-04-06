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

Ggc     = 2
g       = 9.81
rho_w   = 1.0
rho_o   = 0.7
rho_woc = rho_w - rho_o
vc      = 1.0
# mu_o    = 1.0e-3
Kc      = 1.0e-13

# Kc = Ggc * vc * mu_o / (rho_woc * g)
# print(f"{Kc:.2e}")

mu_o = Kc*rho_woc*g/(Ggc*vc)
mu_w = mu_o / 5
print(f"mu_o : {mu_o:.2e}\nmu_w : {mu_w:.2e}")




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