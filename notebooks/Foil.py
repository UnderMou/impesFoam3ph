import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from aux_funcs import *
import scienceplots

def Foil(So, params):
    _,_,_,_,fmoil,floil,epoil = params
    Foil = np.zeros_like(So)
    for i in range(len(So)):
        if So[i] < floil: Foil[i] = 1.0
        elif So[i] > fmoil: Foil[i] = 0.0
        else: Foil[i] = np.power((fmoil-So[i])/(fmoil-floil),epoil)
    
    return Foil

def Fdry(Sw, params):
    _,SF,sfbet,_,_,_,_ = params
    print(SF)
    Fdry = np.zeros_like(Sw)
    for i in range(len(Sw)):
        Fdry = 0.5 + (1.0/np.pi)*np.arctan(sfbet*(Sw-SF))
    return Fdry

if __name__ == '__main__':

    plt.style.use('science')
    plt.rcParams.update({'font.size': 16})
    
    # experiment_name = 'SA_CTscan_FdryFoilFshear_pdr0.2_steady_moreSamples_fixedfmcap_fitted_v2'
    experiment_name = 'SA_CTscan_FdryFoilFshear_pdr0.2_steady_moreSamples_fixedfmcap_fitted_fixedSF'
    experiment_dir = set_experimentDir(experiment_name)
    X_ED = get_samples(experiment_dir)  
    # evidence = [2,9,23]
    evidence = np.arange(0,len(X_ED),1)
    X_ED = pd.DataFrame(X_ED[evidence,:])
    print(X_ED.head())
       

    So = np.linspace(0,1,100)

    Soi = .46
    Foils = []
    for i in range(len(X_ED)):
        Foils.append(Foil(So, X_ED.iloc[i,:].to_numpy()))
    plt.plot(So,Foils[0],c='gray',alpha=0.3,label=f'Samples')
    for i in range(len(X_ED)):
        plt.plot(So,Foils[i],c='gray',alpha=0.3)

    # scatterFoils = []
    # for i in range(len(X_ED)):
    #     scatterFoils.append(Foil([Soi], X_ED.iloc[i,:].to_numpy()))
    # for i in range(len(X_ED)):
    #     plt.scatter(Soi,scatterFoils[i])

    plt.grid()
    plt.legend()
    plt.show()
    plt.close()

    # Fdrys = []
    # for i in range(len(X_ED)):
    #     Fdrys.append(Fdry(So, X_ED.iloc[i,:].to_numpy()))

    # plt.figure(figsize=(5,5))
    # plt.plot(So,Fdrys[0],c='gray',alpha=.3,label=f'Samples')
    # for i in range(len(X_ED)):
    #     plt.plot(So,Fdrys[i],c='gray',alpha=.3)
    
    # plt.xlabel(r"$S_w$")
    # plt.ylabel(r"$F_{dry}$")
    # plt.grid()
    # plt.xlim([0,0.3])
    # plt.legend()
    # plt.show()
    # plt.close()