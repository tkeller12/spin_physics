import numpy as np
from scipy.linalg import expm
from matplotlib.pylab import *
from pyDEER import awg

sigma_x = 0.5*np.r_[
        [
            [0, 1], 
            [1, 0] ]
        ] + 0j

sigma_y = 0.5*np.r_[
        [
            [0,-1j], 
            [1j, 0] ]
        ] + 0j

sigma_z = 0.5*np.r_[
        [
            [1, 0], 
            [0, -1] ]
        ] + 0j

print(sigma_x)
print(sigma_y)
print(sigma_z)


print('commutator test')
print(np.dot(sigma_x,sigma_y) - np.dot(sigma_y,sigma_x))
print(sigma_z * 1j)
print(np.allclose(sigma_z * 1j, np.dot(sigma_x,sigma_y) - np.dot(sigma_y,sigma_x)))

print('Propagator Test')
Px = expm(1j*np.pi/2*sigma_x) # Define Propagator
Py = expm(1j*np.pi/2*sigma_y) # Define Propagator


sigma = sigma_z
print('Px')
print(Px)
print('Py')
print(Py)

#out = P @ sigma_z @ P.conj()

print('Sigma y')
#print(sigma_y)
print('output')
#print(out)




