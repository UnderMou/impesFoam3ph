import pandas as pd
import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import scienceplots

# def sigma_from_DP(DP):
#     """
#     Compute lognormal standard deviation (sigma)
#     from a given DP coefficient.
#     """

#     def func(sigma):
#         return 1.0 - np.exp(-1.5 * sigma**2 - sigma) - DP

#     sigma_initial = 0.5
#     sigma_solution = fsolve(func, sigma_initial)[0]

#     return sigma_solution


def permeability_from_DP(DP, k_mean, n=100000):
    """
    Generate permeability distribution using
    DP definition from log-normal probability plot.
    """

    # standard deviation in log10 space
    sigma_logk = -np.log10(1.0 - DP)

    # mean in log space
    mu_logk = np.log10(k_mean)

    # sample
    logk = np.random.normal(mu_logk, sigma_logk, n)

    return 10**logk, 10**mu_logk, 10**sigma_logk

def sigma_from_DP(DP):
    """
    Compute lognormal standard deviation (sigma)
    from a given DP coefficient.
    """

    def func(sigma):
        return 1.0 - np.exp(-1.5 * sigma**2 - sigma) - DP
        # return 1.0 - np.exp(-0.5 * sigma**2 - sigma) - DP

    sigma_initial = 0.5
    sigma_solution = fsolve(func, sigma_initial)[0]

    return sigma_solution


def permeability_distribution(DP, k_mean, n_samples=10000):
    """
    Generate permeability distribution from DP coefficient.
    """

    sigma = sigma_from_DP(DP)

    # compute mu from mean permeability
    mu = np.log(k_mean) - 0.5 * sigma**2

    # generate permeability field
    k = np.random.lognormal(mean=mu, sigma=sigma, size=n_samples)

    return k, mu, sigma

if __name__ == "__main__":

    plt.style.use('science')

    plt.rcParams["font.size"] =  22
    
    np.random.seed(1)

    DP = [0.05, 0.4, 0.8, 0.95]
    k_mean = 100.0  # mD
    
    # print(f"mean k = {k_mean*9.869233e-16} m2")

    # outputs = [permeability_from_DP(dp, k_mean) for dp in DP]
    outputs = [permeability_distribution(dp, k_mean) for dp in DP]
    k, mu, sigma = zip(*outputs)

    # for i,_ in enumerate(k):
    #     k[i][:] *= 9.869233e-16
        

    fig, ax = plt.subplots(1, len(DP), figsize=(20,6), sharey=False)

    # global limits (important!)
    k_all = np.concatenate(k)
    xlims = [k_all.min(), k_all.max()]

    # log-spaced bins
    logbins = np.logspace(
        np.log10(xlims[0]),
        np.log10(xlims[1]),
        80
    )

    for i in range(len(DP)):

        ax[i].hist(k[i], bins=logbins)
        ax[i].set_xscale("log")

        ax[i].set_xlabel("$k$ [mD]")
        ax[i].set_title(f"DP = {DP[i]}")
        ax[i].grid(True, which="both")

        ax[i].set_xlim(xlims)

    ax[0].set_ylabel("frequency")

    plt.tight_layout()
    plt.savefig('DP_distribution.pdf', dpi=300)
    plt.show()
    plt.close()

    fig, ax = plt.subplots(1, len(DP), figsize=(20,6), sharey=False)


    k_all = np.concatenate(k)

    xlims = [k_all.min(), k_all.max()]

    logbins = np.logspace(
        np.log10(xlims[0]),
        np.log10(xlims[1]),
        80
    )

    for i in range(len(DP)):

        ax[i].hist(np.log10(k[i]), bins=60)

        ax[i].set_xlabel(r"$\log_{10}(k)$")
        ax[i].set_title(f"DP = {DP[i]}")
        ax[i].grid(True)
    
    ax[0].set_ylabel("frequency")

    plt.tight_layout()
    plt.show()