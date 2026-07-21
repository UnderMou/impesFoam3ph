import numpy as np

def krw(S):
    return 0.302 * S**2

def krg(S):
    return 0.004 * (1-S)**2

def FM(Sw,Cs):
    fmmob = 293.27
    epdry = 359.33
    fmdry = 0.437
    fmsurf = 2.0
    epsurf = 1.0

    F2 = 1/2 + np.arctan(epdry*(Sw - fmdry))/np.pi

    if Cs >= fmsurf: F1 = 1.0
    else: F1 = np.power(Cs/fmsurf, epsurf)

    return 1/(1 + fmmob*F1*F2)

if __name__ == '__main__':

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

    # Fractional flow - INPUT (L state)
    muw = 1e-3
    mug = 5e-5

    krg_f = krg(Sw_L_norm) * FM(Sw_L, C_L)

    fw = krw(Sw_L_norm) / (krw(Sw_L_norm) + muw/mug*krg_f)
    print('fw = ', fw, ' | fg = ', 1-fw)