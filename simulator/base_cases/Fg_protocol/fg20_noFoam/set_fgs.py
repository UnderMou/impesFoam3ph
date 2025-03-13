import numpy as np

U = 5.84748220e-06 + 2.33899288e-05
fgs = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])  
fws = 1.0 - fgs

Uas = fgs * U
Ubs = fws * U

for i in range(len(fgs)):
    print('fg = ', fgs[i], ' | Ua = ', Uas[i], ' | Ub = ', Ubs[i])