# Parameters
Swc = 0.437
Sgr = 0.293
phi = 0.21
L = 1.0
Cmax = 2

# Normalized variables
C_L_norm = 0.5
C_R_norm = 0
C_L = C_L_norm * Cmax
C_R = C_R_norm * Cmax
print("C_L = ", C_L, '|', 'C_R = ', C_R)

Sw_L_norm = 0.99
Sw_R_norm = 0.01
Sw_L = Sw_L_norm * (1-Swc-Sgr) + Swc
Sw_R = Sw_R_norm * (1-Swc-Sgr) + Swc
print("Sw_L = ", Sw_L, '|', 'Sw_R = ', Sw_R)


# Control
tD = 0.3
u = 1.0
tf = tD*(1-Swc-Sgr)*phi*L / u
print('tf = ', tf,'[s]')