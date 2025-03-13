import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

def calc_Se(params,Sb):
    Sb_min = params['Sb_min']
    Sa_min = params['Sa_min']

    return np.divide(Sb - Sb_min, 1-Sa_min - Sb_min)

def calc_krb(params, Sb):
    krb_max = params['krb_max']
    b = params['b']

    Se = calc_Se(params,Sb)

    return krb_max * np.power(Se, b)

def calc_kra(params, Sb):
    kra_max = params['kra_max']
    a = params['a']

    Se = calc_Se(params,Sb)

    return kra_max * np.power(1.0-Se, a)

def calc_Fb(params, Sb):
    krb = calc_krb(params, Sb)
    kra = calc_kra(params, Sb)

    mu_b = params['mu_b']
    mu_a = params['mu_a']

    Fb = (krb/mu_b) / (krb/mu_b + kra/mu_a)
    return Fb

def func(x, params):
    Fb_inj = params['Fb_inj']

    mob_b = calc_krb(params, x) / params['mu_b']
    mob_a = calc_kra(params, x) / params['mu_a']

    return Fb_inj - mob_b / (mob_b+mob_a)

if __name__ == "__main__":

    params = {
        'Fb_inj':0.2,
        'U_in':1.0e-5,
        'Sb_min':0.1,
        'Sa_min':0.15,
        'krb_max':1.0,
        'a':8,
        'kra_max':0.5,
        'b':3,
        'mu_b': 0.001,
        'mu_a': 1.0e-5,
        'K' : 1.0e-13
    }

    Sb_in = optimize.newton(func, 0.5, args=(params,))
    print("Sb_in: ", Sb_in)
    print("Fb_in: ", params['Fb_inj'], " | Fb(Sb_in): ", calc_Fb(params,Sb_in))

    # inlet
    U_in = params['U_in']
    K = params['K']
    alpha = K * (calc_krb(params,Sb_in)/params['mu_b'] + calc_kra(params,Sb_in)/params['mu_a'])
    gradP = U_in / alpha
    print('gradP = ', gradP)

    Ub_in = U_in * params['Fb_inj']
    Ua_in = U_in - Ub_in
    print("Ub_in: ", Ub_in,"\nUa_in: ", Ua_in, '\nU: ', Ub_in+Ua_in)

    exit()
    
    # plot fw
    Sw = np.linspace(params['Swc'], 1.0-params['Sor'], 100, endpoint=False)
    fw = calc_fw(params, Sw)
    plt.plot(Sw, fw, label='fw_python')
    SwOF = np.array([0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.6,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.7,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.8,0.81,0.82,0.83,0.84,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85])
    fwOF = np.array([0,0,0,0,0,0,0,0,0,0,0,5.19958e-15,1.38655e-12,3.7037e-11,3.85806e-10,2.39955e-09,1.07727e-08,3.86297e-08,1.17533e-07,3.1548e-07,7.67226e-07,1.72292e-06,3.62318e-06,7.21161e-06,1.36991e-05,2.49994e-05,4.40604e-05,7.53249e-05,0.000125362,0.000203723,0.000324097,0.000505844,0.000776033,0.00117211,0.00174536,0.00256542,0.00372598,0.00535204,0.00760882,0.0107126,0.0149437,0.0206601,0.0283129,0.0384588,0.0517676,0.069019,0.0910786,0.118845,0.153155,0.194651,0.243606,0.299746,0.362121,0.429079,0.498388,0.567511,0.633959,0.695628,0.751027,0.799355,0.840446,0.87462,0.902513,0.924924,0.942692,0.956622,0.967438,0.975762,0.982118,0.986933,0.990553,0.993253,0.995247,0.996706,0.997761,0.998513,0.99904,0.999402,0.999644,0.999801,0.999897,0.999953,0.999982,0.999995,0.999999,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
    plt.scatter(SwOF, fwOF, label='fw_OF')
    plt.legend()
    plt.show()

    # plot rel.perm.
    krw = calc_krw(params, Sw)
    kro = calc_kro(params, Sw)
    plt.plot(Sw, krw, label='krw_python')
    plt.plot(Sw, kro, label='kro_python')
    plt.legend()
    plt.show()