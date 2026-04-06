import uqtopus as uqt
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import ternary
import matplotlib
from scipy import integrate

dir_path = os.path.dirname(os.path.realpath(__file__))


data_single = uqt.parse_openfoam_case(
    case_dir=dir_path,
    variables=["Sa", "Sb", "Sc", "Fa", "Fb", "p"],
    time_dirs=[f"{i:g}" for i in np.arange(25, 5000, 25)]
)
Sc = data_single.Sc[-1,:]

Soi = 1 - 0.0131 - 0.527

dPs = []
Fcs = []
GOR = []
I = []
t = []
for j in range(len(data_single.time)):

    Sc = data_single.Sc[j,:]
    x = np.linspace(0,0.4,len(Sc))

    I1 = integrate.simpson(Sc,x=x)
    print('t = ', np.asarray(data_single.time[j]), ' | Int(Sc) = ', I1)

    I.append(I1)
    t.append(np.asarray(data_single.time[j]))
    
    Fa = np.asarray(data_single.Fa[j,:])[-1]
    Fb = np.asarray(data_single.Fb[j,:])[-1]
    Fc = 1 - Fa- Fb
    GOR.append(Fa/Fc)

    Fcs.append(Fc)
    
    p_in = np.asarray(data_single.p[j,:])[0]
    p_out = np.asarray(data_single.p[j,:])[-1]
    dP = p_in - p_out
    dPs.append(dP)
    



plt.plot(t, I, c='b',alpha=0.5)
plt.title('Integral of Sc over Time')
plt.xlabel('Time')
plt.ylabel('Integral of Sc')
plt.grid(True)
plt.legend(['Integral of Sc'])
plt.show()
plt.close()

RF = 1 - (np.array(I)/(Soi*0.4))
plt.plot(t, RF*100, c='b',alpha=0.5)
plt.title('Recovery Factor over Time')
plt.xlabel('Time')
plt.ylabel('Recovery Factor (%)')
plt.grid(True)
plt.legend(['Recovery Factor'])
plt.show()
plt.close()

plt.plot(t, GOR, c='b',alpha=0.5)
plt.title('Gas-Oil Ratio (GOR) over Time')
plt.xlabel('Time')
plt.ylabel('Gas-Oil Ratio (GOR)')
plt.grid(True)
plt.legend(['Gas-Oil Ratio'])
plt.tight_layout()  # Adjust layout to make room for labels
plt.show()
plt.close()

plt.plot(t, Fcs, c='b', alpha=0.5)
plt.title('Oil fractional flow')
plt.xlabel('Time')
plt.ylabel('fo')
plt.grid(True)
plt.legend(['fo'])
plt.tight_layout()
plt.show()

plt.plot(t, dPs, c='b', alpha=0.5)
plt.title('Pressure Drop (dP) over Time')
plt.xlabel('Time')
plt.ylabel('Pressure Drop (MPa)')
plt.grid(True)
plt.legend(['Pressure Drop'])
plt.tight_layout()
plt.show()
