import math

Sb_min = 0.1
Sa_min = 0.15
kra_max = 0.5
a = 8.0
krb_max = 1.0
b = 3.0

Sb = 0.4375
Se = (Sb - Sb_min) / (1-Sa_min-Sb_min)

krb = krb_max * pow(Se, b)
kra = kra_max * pow(1-Se, a) 

print(Sb, kra, krb)