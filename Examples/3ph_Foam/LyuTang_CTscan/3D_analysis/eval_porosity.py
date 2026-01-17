import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Define the normal distribution function
def normal_dist(x, mu, sigma, amplitude):
    return amplitude * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

if __name__ == '__main__':
    # Load the data
    data = pd.read_csv('porosity_freqHist.csv')
    x_data = data.iloc[:, 0].to_numpy()
    y = data.iloc[:, 1].to_numpy()

    # Fit the normal distribution to the data
    popt, pcov = curve_fit(normal_dist, x_data, y, p0=[np.mean(x_data), np.std(x_data), max(y)])

    # Extract the fitted parameters
    mu, sigma, amplitude = popt

    # Generate the fitted y values
    y_fitted = normal_dist(x_data, mu, sigma, amplitude)

    # Plot the original data and the fitted distribution
    plt.figure(figsize=(8, 6))
    plt.scatter(x_data, y, label='Original Data', color='blue')
    plt.plot(x_data, y_fitted, label=f'Fitted Normal Dist\n$\\mu={mu:.2f}, \\sigma={sigma:.2f}$', color='red', lw=2)
    plt.xlabel('x (Porosity)')
    plt.ylabel('Frequency')
    plt.title('Normal Distribution Fit')
    plt.legend()
    plt.grid()
    plt.show()

    # Print the fitted parameters
    print(f"Fitted Parameters:\nMean (μ): {mu:.2f}\nStandard Deviation (σ): {sigma:.2f}\nAmplitude: {amplitude:.2f}")


    # Number of samples to generate
    num_samples = 38400

    # Generate samples from the normal distribution using the fitted mean and std dev
    samples = np.random.normal(loc=mu, scale=sigma, size=num_samples)

    # Plot the histogram of the samples for verification
    plt.figure(figsize=(8, 6))
    plt.hist(samples, bins=30, density=True, alpha=0.6, color='g', label='Sampled Data Histogram')

    # Overlay the fitted normal PDF for comparison
    x_fine = np.linspace(min(samples), max(samples), 500)
    pdf_fitted = normal_dist(x_fine, mu, sigma, amplitude)  # Use the fitted function
    plt.plot(x_fine, pdf_fitted, 
            'r-', label='Fitted Normal PDF')

    # Add labels and legend
    plt.xlabel('x (Porosity)')
    plt.ylabel('Density')
    plt.title('Sampled Data vs. Fitted PDF')
    plt.legend()
    plt.grid()
    plt.show()

    # Write a porosity file from samples
    np.savetxt('porosity_samples.csv', samples, delimiter=',', fmt='%f')

    # Write a permeability file from fitted PDF
    mu_K = 2.783123706e-12
    sigma_K = (sigma/mu) * mu_K
    samples = np.random.normal(loc=mu_K, scale=sigma_K, size=num_samples)
    np.savetxt('permeability_samples.csv', samples, delimiter=',', fmt='%.6e')
    plt.figure(figsize=(8, 6))
    plt.hist(samples, bins=30, density=True, alpha=0.6, color='g', label='Sampled Data Histogram')
    plt.xlabel('x (Permeability)')
    plt.ylabel('Density')
    plt.legend()
    plt.grid()
    plt.show()
