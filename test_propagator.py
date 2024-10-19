import numpy as np
from scipy.linalg import expm
from matplotlib.pylab import *

# Define spin operators
Sx = 0.5*np.r_[
        [
            [0, 1], 
            [1, 0] ]
        ] + 0j

Sy = 0.5*np.r_[
        [
            [0,-1j], 
            [1j, 0] ]
        ] + 0j

Sz = 0.5*np.r_[
        [
            [1, 0], 
            [0, -1] ]
        ] + 0j

print(Sx)
print(Sy)
print(Sz)


print('commutator test')
print(np.dot(Sx,Sy) - np.dot(Sy,Sx))
print(Sz * 1j)
print(np.allclose(Sz * 1j, np.dot(Sx,Sy) - np.dot(Sy,Sx)))

print('Propagator Test')
Px = expm(1j*np.pi/2*Sx) # Define Propagator
Py = expm(1j*np.pi/2*Sy) # Define Propagator


sigma = Sz
print('Px')
print(Px)
print('Py')
print(Py)

#out = Px @ Sz @ Px.T.conj()
out = Py @ Sz @ Py.T.conj()

print('Sigma y')
print(Sy)
print('output')
print(out)




