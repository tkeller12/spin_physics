import numpy as np
from scipy.linalg import expm
from matplotlib.pylab import *

sigma_x = 0.5*np.r_[[[0, 1],[1, 0]]]
sigma_y = 0.5*np.r_[[[0,-1j],[1j, 0]]]
sigma_z = 0.5*np.r_[[[1, 0],[0, -1]]]

print(sigma_x)
print(sigma_y)
print(sigma_z)


print('commutator test')
print(np.dot(sigma_x,sigma_y) - np.dot(sigma_y,sigma_x))
print(sigma_z * 1j)

omega = 10.
dt = 0.1

coil = sigma_x + 1j*sigma_y

H = omega * sigma_z

P = expm(1j*H*dt)

sigma = sigma_x
M_list = []

for ix in range(100):
    sigma = np.dot(np.dot(P,sigma),np.conjugate(P))
    M = np.trace(np.dot(coil,sigma))
    M_list.append(M)


figure()
plot(M_list)
show()
