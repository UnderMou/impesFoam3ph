import pandas as pd
import matplotlib.pyplot as plt
import scienceplots
import numpy as np

plt.style.use('science')

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['font.size'] = 16

data = pd.read_csv("production.txt", delimiter=',')
print(data.head())
print(data.tail())

t = data['TimeStep'].to_numpy()[1:]
qo = data['qo'].to_numpy()[1:]
qw = data['qw'].to_numpy()[1:]

Qw = np.zeros_like(qw)
Qo = np.zeros_like(qo)

for i in range(len(Qw)):
    Qw[i] = np.sum(qw[:i])
    Qo[i] = np.sum(qo[:i])

# plt.figure(figsize=(7,7))
# plt.plot(t,qw,c='b',label=r"$Q_w$")
# plt.plot(t,qo,c='k',label=r"$Q_o$")
# plt.grid()
# plt.legend()
# # plt.xlim([0,t[-1]+3000])
# # plt.xlim([0,300000])
# # plt.ylim([0,4])
# plt.xlabel('time [s]')
# plt.ylabel('production []')
# plt.show()



fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15,5))

ax[0].plot(t,qw,c='b',label=r"$q_w$")
ax[0].plot(t,qo,c='k',label=r"$q_o$")
ax[0].grid()
ax[0].legend()
ax[0].set_xlabel('time [s]')
ax[0].set_ylabel(r'production [$m^3/s$]')

ax[1].plot(t,Qw,c='b',label=r"$Q_w$")
ax[1].plot(t,Qo,c='k',label=r"$Q_o$")
ax[1].grid()
ax[1].legend()
ax[1].set_xlabel('time [s]')
ax[1].set_ylabel(r'cumulative production [$m^3$]')

plt.tight_layout()
plt.savefig("Production.pdf", dpi=300)
plt.show()