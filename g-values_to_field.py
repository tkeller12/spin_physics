import numpy as np
from matplotlib.pylab import *
import dnplab as dnp

# typical g-values for MTSSL:
# gx = 2.0083-2.0091, gy = 2.0061, gz = 2.0022
g = np.r_[2.0083, 2.0061, 2.0022]

mu_B = dnp.physical_constants['Bohr magneton'][0]
hbar = dnp.hbar

f = 34.4e9

gamma = g * mu_B / hbar

B = 2*np.pi * f / gamma

print(gamma)
print(B)
