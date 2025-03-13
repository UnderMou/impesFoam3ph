import numpy as np
import matplotlib.pyplot as plt

Sw = np.linspace(0,1,100)
Sg = Sw

pcow = 5*(2-Sw)*(1-Sw)
pcgo = (2-Sg)*(1-Sg)

fig, axs = plt.subplots(1, 2, figsize=(10,5))

axs[0].plot(Sw,np.gradient(pcow)*np.gradient(Sw))
axs[1].plot(Sg,np.gradient(pcgo)*np.gradient(Sg))

plt.show()