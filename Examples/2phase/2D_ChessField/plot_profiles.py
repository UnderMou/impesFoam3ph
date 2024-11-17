import pandas as pd
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')
# plt.rcParams.update({'font.size': 16})

# profiles = ['CUBICK','UPWIND','limitedLINEAR','SuperBee','TOPUS']
profiles = ['UPWIND','TOPUS','limitedLINEAR','CUBICK','SuperBee']

plt.figure(figsize=(5,4))
path = './profiles/'
for i in range(len(profiles)):
    data = pd.read_csv(path+profiles[i]+'.csv')
    x = data['Points:0']
    y = data['Sb']
    plt.plot(x,y,label=profiles[i])
plt.legend()
plt.xlabel(r'$x$ [m]')
plt.ylabel(r'$S_g$ [-]')
plt.grid()
plt.savefig('schemes_comparative.pdf',dpi=300)
plt.show()