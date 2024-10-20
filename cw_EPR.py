# Timothy Keller
# S = 1/2, I = 1/2
# Spin 1/2 electron coupled to spin 1/2 nuclei

import numpy as np
from scipy.linalg import expm
from matplotlib.pylab import *
from matplotlib import cm

sigma_x = 0.5*np.r_[[[0, 1],[1, 0]]]
sigma_y = 0.5*np.r_[[[0,-1j],[1j, 0]]]
sigma_z = 0.5*np.r_[[[1, 0],[0, -1]]]
Identity = np.eye(2)


Sx = np.kron(sigma_x, Identity)
Sy = np.kron(sigma_y, Identity)
Sz = np.kron(sigma_z, Identity)

Ix = np.kron(Identity, sigma_x)
Iy = np.kron(Identity, sigma_y)
Iz = np.kron(Identity, sigma_z)

SxIx = np.kron(sigma_x,sigma_z)

SxIx2 = np.dot(Sx,Iz)

coil = Sx + 1j*Sy

print(SxIx)
print(SxIx2)
print(np.allclose(SxIx,SxIx2))
omega_S = 1.76e11 # rad / (s * T)
omega_I = 267.522e6 # rad / (s * T)
Aiso = 2*np.pi * 50.e6 # Isotropic Hyperfine coupling rad / s

B0 = 0.35# T

H = omega_S/(2.*np.pi)*B0*Sz + omega_I/(2.*np.pi)*B0*Iz + Aiso * np.dot(Sz,Iz)
#H = omega_S/(2.*np.pi)*B0*Sz + omega_I/(2.*np.pi)*B0*Iz + Aiso * (np.dot(Sx,Ix) + np.dot(Sy,Iy) + np.dot(Sz,Iz))

print('Hamiltonian:')
print(H)
out = np.linalg.eig(H)

E = out[0]
E_vectors = out[1]

vect1 = E_vectors[:,0]
vect2 = E_vectors[:,3].T

V = np.dot(np.dot(vect1,coil),vect2)

print('V')
print(V)




print('Eigenvalues')
print(E)
print('Eigenvectors')
print(E_vectors)

E12 = E[0] - E[1]
E34 = E[2] - E[3]
E13 = E[0] - E[2]
E24 = E[1] - E[3]
print(E12)
print(E34)
print(E13)
print(E24)

print('Nuclear')
print('%0.05f MHz'%(E12 / 1e6))
print('%0.05f MHz'%(E34 / 1e6))
print('Electron')
print('%0.05f GHz'%(E13 / 1e9))
print('%0.05f GHz'%(E24 / 1e9))

matshow(abs(H), cmap = cm.jet)
title('Hamiltonian')
show()
