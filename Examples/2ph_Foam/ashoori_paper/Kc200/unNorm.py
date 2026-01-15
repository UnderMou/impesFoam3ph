import numpy as np

def krw(Sw,Swc,Sgr):
    S = (Sw-Swc)/(1-Swc-Sgr)
    return 0.2 * S**4.2

def krg(Sw,Swc,Sgr):
    S = (1-Sw-Sgr)/(1-Swc-Sgr)
    return 0.94 * S**1.3

def MRF(Sw,nD):
    Cmrf = 18500

    MRF = Cmrf*nD + 1   

    return MRF

if __name__ == '__main__':

    # Parameters
    Swc = 0.2
    Sgr = 0.18

    # Fractional flow - INPUT (L state)
    u_inj = 2.93e-8
    Sw_inj = 0.372
    nD_inj = 0.664
    muw = 1.0e-3
    mug = 2.0e-5

    krg_f = krg(Sw_inj,Swc,Sgr) / MRF(Sw_inj,nD_inj)

    fw = krw(Sw_inj,Swc,Sgr) / (krw(Sw_inj,Swc,Sgr) + muw/mug*krg_f)
    print('fw = ', fw, ' | fg = ', 1-fw)
    print('uw = ', u_inj*fw, ' | ug = ', u_inj*(1-fw))