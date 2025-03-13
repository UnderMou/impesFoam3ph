import numpy as np
import math

p = 8
g = 4
alpha = 3
trainingSet_size = 0.75

Nsamples = alpha * math.factorial(p + g) / (math.factorial(p) * math.factorial(g))
Nsamples = Nsamples / trainingSet_size
print("Nsamples: ", int(Nsamples))