import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv('perm_inlet.csv')
# print(data.head())
y = data['Points:1'].to_numpy()
K = data['K'].to_numpy()

plt.plot(y,K)
plt.show()

print(K.mean(), np.median(K))