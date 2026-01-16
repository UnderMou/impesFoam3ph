import numpy as np

def MRF(params, Sa, Sb):

    fmmob = params['fmmob']
    SF = params['SF']
    sfbet = params['sfbet']
    fmoil = params['fmoil']
    floil = params['floil']
    epoil = params['epoil']

    # Fdry
    Fdry = 0.5 + (1.0 / np.pi) * np.arctan(sfbet * (Sb - SF))
    
    # Foil
    Sc = 1.0 - Sa - Sb
    Foil = np.zeros_like(Sc)
    
    if fmoil <= Sc and Sc <= 1.0 - params['Sa_min'] - params['Sb_min']:
        Foil = 0.0
    elif floil <= Sc and Sc < fmoil:
        Foil = np.power((fmoil-Sc)/(fmoil-floil),epoil)
    else:
        Foil = 1.0

    MRF = 1.0 / (1.0 + fmmob*Fdry*Foil)
    # print(MRF)

    return MRF

def calc_kra(params, Sa, Sb):

    Se_a = (Sa - params['Sa_min']) / ( 1 - params['Sa_min'] - params['Sb_min'] - params['Sc_min'] )

    return params['kra_max'] * np.power(Se_a, params['a']) * MRF(params, Sa, Sb)

def calc_krb(params, Sb):

    Se_b = (Sb - params['Sb_min']) / ( 1 - params['Sa_min'] - params['Sb_min'] - params['Sc_min'] )

    return params['krb_max'] * np.power(Se_b, params['b'])

def calc_krc(params, Sa, Sb):

    Se_a = (Sa - params['Sa_min']) / ( 1 - params['Sa_min'] - params['Sb_min'] - params['Sc_min'] )
    Se_b = (Sb - params['Sb_min']) / ( 1 - params['Sa_min'] - params['Sb_min'] - params['Sc_min'] )

    return params['krc_max'] * np.power(1 - Se_a - Se_b, params['c'])

def calc_Fa(params, Sa, Sb):
    kra = calc_kra(params, Sa, Sb) 
    krb = calc_krb(params, Sb)
    krc = calc_krc(params, Sa, Sb)
    
    mu_c = params['mu_c']
    mu_b = params['mu_b']
    mu_a = params['mu_a']

    return (kra/mu_a) / (krb/mu_b + kra/mu_a + krc/mu_c)

def calc_Fb(params, Sa, Sb):
    kra = calc_kra(params, Sa, Sb) 
    krb = calc_krb(params, Sb)
    krc = calc_krc(params, Sa, Sb)
    
    mu_c = params['mu_c']
    mu_b = params['mu_b']
    mu_a = params['mu_a']

    return (krb/mu_b) / (krb/mu_b + kra/mu_a + krc/mu_c)

def calc_Fc(params, Sa, Sb):
    kra = calc_kra(params, Sa, Sb) 
    krb = calc_krb(params, Sb)
    krc = calc_krc(params, Sa, Sb)
    
    mu_c = params['mu_c']
    mu_b = params['mu_b']
    mu_a = params['mu_a']

    return (krc/mu_c) / (krb/mu_b + kra/mu_a + krc/mu_c)

if __name__ == "__main__":

    params = {
        'Sa_min':0.0,
        'Sb_min':0.1,
        'Sc_min':0.1,
        'kra_max':1.0,
        'krb_max':1.0,
        'krc_max':1.0,
        'a':2,
        'b':2,
        'c':2,
        'mu_a': 1.0e-5,
        'mu_b': 1.0e-3,
        'mu_c': 5.0e-3,
        'fmmob' : 2000.0,
        'SF' : 0.3,
        'sfbet' : 32000.0,
        'fmoil' : 0.3,
        'floil' : 0.1,
        'epoil' : 3.0
    }

    # BC definition
    U_in = 1.4111e-5
    Sa_in = 0.64
    Sb_in = 0.26

    fa = calc_Fa(params, Sa_in, Sb_in)
    fb = calc_Fb(params, Sa_in, Sb_in)
    fc = calc_Fc(params, Sa_in, Sb_in)

    Ua = fa * U_in
    Ub = fb * U_in
    Uc = fc * U_in

    print('U : ', U_in)
    print('fa : ', fa, ' | fb : ', fb, ' | fc : ', fc)
    print('Ua : ', Ua, ' | Ub : ', Ub, ' | Uc : ', Uc)