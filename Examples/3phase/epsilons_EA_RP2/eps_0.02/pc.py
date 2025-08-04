import numpy as np
import matplotlib.pyplot as plt

def pcow(sw,epsilon):
    return -5*epsilon*(2-sw)*(1-sw)

def pcgo(sg,epsilon):
    return 5*epsilon*(2-sg)*(1-sg)

epsilon = 0.01
sw = np.linspace(0,1,1000)
sg = sw
pcow_v = pcow(sw,epsilon)
pcgo_v = pcgo(sg,epsilon)

fig, ax = plt.subplots(1,2,figsize=(7,3))

ax[0].plot(sw,pcow_v)
ax[1].plot(sg,pcgo_v)

plt.tight_layout()
plt.show()

print(np.max(pcow_v)-np.min(pcow_v), np.mean(pcow_v))