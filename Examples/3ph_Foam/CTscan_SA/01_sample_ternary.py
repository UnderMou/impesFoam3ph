import numpy as np
import matplotlib.pyplot as plt
from uqpylab import sessions, display_util
import seaborn as sns
import pandas as pd
from scipy.special import comb

def Foil (So, fmoil, floil, epoil):
    
    F = np.zeros((len(fmoil), len(So)))

    for i in range(len(fmoil)):
        fmoil_i = fmoil[i]
        floil_i = floil[i]
        epoil_i = epoil[i]

        for j in range(len(So)):
            if So[j] < floil_i: F[i, j] = 1.0
            elif So[j] > fmoil_i: F[i, j] = 0.0
            else: F[i, j] = ((fmoil_i - So[j]) / (fmoil_i - floil_i)) ** epoil_i

    return F


if __name__ == "__main__":

    # User defined ===========
    case_name = "CTscanSA"
    # ========================

    path2save = f"casesPrep/{case_name}/"
    
    ############ UQ[py]Lab - INPUT #############

    # Initialize common plotting parameters
    display_util.load_plt_defaults()
    uq_colors = display_util.get_uq_color_order()

    # Start the session
    myToken = '270fe799b4e4e4f004da79cec3dffb07cf9da78a' # The user's token to access the UQCloud API
    UQCloud_instance = 'https://uqcloud.ethz.ch' # The UQCloud instance to use
    mySession = sessions.cloud(host=UQCloud_instance, token=myToken)
    # (Optional) Get a convenient handle to the command line interface
    uq = mySession.cli
    # Reset the session
    mySession.reset()

    # Set random seed for reproducibility
    uq.rng(0,'twister')

    # Reference values
    fmmobR = 50000.0
    sfdryR = 0.215
    sfbetR = 19950.0
    epcapR = 1.321
    fmoilR = 0.823
    floilR = 0.295
    epoilR = 3.827

    # # Uncertainty characterization - COV de outros papers
    # cov_fmmob = 0.09066994                  # WRR, 2025 
    # cov_SF = 0.001030595                    # WRR, 2025 
    # cov_sfbet = 0.274304185                 # WRR, 2025
    # cov_epcap = 0.096058532                 # WRR, 2025
    # cov_fmoil = 0.298                       # deMiranda, 2025
    # cov_floil = 0.271                       # deMiranda, 2025
    # cov_epoil = 0.358                       # deMiranda, 2025

    # Uncertainty characterization - COV de 10%
    cov_fmmob = 0.1 
    cov_SF = 0.1 
    cov_sfbet = 0.1
    cov_epcap = 0.1
    cov_fmoil = 0.1
    cov_floil = 0.1
    cov_epoil = 0.1

    dists = {
    'num_vars': 7,
    'names': ['fmmob', 'sfdry', 'sfbet', 'epcap', 'fmoil', 'floil', 'epoil'],
    'bounds': [[fmmobR*(1-cov_fmmob/2), fmmobR*(1+cov_fmmob/2)],
               [sfdryR*(1-cov_SF/2), sfdryR*(1+cov_SF/2)],
               [sfbetR*(1-cov_sfbet/2), sfbetR*(1+cov_sfbet/2)],
               [epcapR*(1-cov_epcap/2), epcapR*(1+cov_epcap/2)],
               [fmoilR*(1-cov_fmoil/2), fmoilR*(1+cov_fmoil/2)],
               [floilR*(1-cov_floil/2), floilR*(1+cov_floil/2)],
               [epoilR*(1-cov_epoil/2), epoilR*(1+cov_epoil/2)],
              ]
    }

    # Set marginals
    InputOpts = {
        "Marginals": [
            {"Name": dists['names'][0],
            "Type": "Uniform",
            "Parameters": dists['bounds'][0]
            },
            {"Name": dists['names'][1],
            "Type": "Uniform",
            "Parameters": dists['bounds'][1]
            },
            {"Name": dists['names'][2],
            "Type": "Uniform",
            "Parameters": dists['bounds'][2]
            },
            {"Name": dists['names'][3],
            "Type": "Uniform",
            "Parameters": dists['bounds'][3]
            },
            {"Name": dists['names'][4],
            "Type": "Uniform",
            "Parameters": dists['bounds'][4]
            },
            {"Name": dists['names'][5],
            "Type": "Uniform",
            "Parameters": dists['bounds'][5]
            },
            {"Name": dists['names'][6],
            "Type": "Uniform",
            "Parameters": dists['bounds'][6]
            }
        ]
    }

    myInput = uq.createInput(InputOpts)

    uq.print(myInput)
    # uq.display(myInput)

    # Number of samples for PCE/Sobol
    d = 7       # parameters
    p = 3       # degree
    n_terms = int(comb(d + p, p))

    alpha = 5
    nSamples = alpha*n_terms
    print(nSamples)

    X_ED = uq.getSample(myInput, max(nSamples,1000))

    So=np.linspace(0,1,100)
    Foil_vec = Foil(
        So,
        fmoil=X_ED[:,4],
        floil=X_ED[:,5],
        epoil=X_ED[:,6]
    )
    for i in range(Foil_vec.shape[0]):
        plt.plot(So, Foil_vec[i,:], c='gray', alpha=0.4)
    plt.savefig(f"{path2save}Foil.pdf", dpi=300)

    # plot
    df = pd.DataFrame(
        X_ED,
        columns=dists['names']
    )
    print(df.head())
    print(df.shape)

    floilBigger = df["floil"] > df["fmoil"]
    print('Any floil > fmoil: ', any(floilBigger))
    print(df.describe())

    # Pairplot with density-style histograms
    sns.pairplot(
        df,
        diag_kind='kde',      # diagonal density curves
        corner=True,          # lower triangle only
        plot_kws={
            's': 12,
            'alpha': 0.5
        }
    )

    plt.tight_layout()
    plt.savefig(f"{path2save}input_samples.pdf", dpi=300)

    df.to_csv(f"{path2save}input_samples.csv", index=False)