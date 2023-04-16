import numpy as np
from scipy.linalg import expm
from matplotlib.pylab import *

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

omega = 100000. # frequency offset from carrier, Hz
#tp = 0.0001 # Pulse Length, s
B1 = 20000
pts = 1024 # Points in FID

tp = np.pi/2 / B1

omega_array = np.r_[-omega/2:omega/2:1j*pts]

coil = sigma_x + 1j*sigma_y # Detection Operator (NMR Coil)

#H = 2*np.pi * omega * sigma_z # Calculate Hamiltonian (only Zeeman)

#P = expm(1j*H*dt) # Define Propagator

#sigma = sigma_z # Initial Density Matrix (After 90-pulse)

M_list = []

for ix in range(pts):
    sigma = sigma_z # Initial Density Matrix
    # re-calculate spin hamiltonian for offset
#    H = 2*np.pi * omega_array[ix] * sigma_z + np.pi/2 * sigma_y # Calculate Hamiltonian (only Zeeman)

    H = 2*np.pi * tp * omega_array[ix] * sigma_z + B1 * tp * sigma_x # Calculate Hamiltonian (only Zeeman)
    P = expm(1j*H) # Define Propagator
#    sigma = np.dot(np.dot(P,sigma),np.conjugate(P)) # Propagate Density Matrix
    sigma = np.dot(np.dot(P,sigma),P.conj()) # Propagate Density Matrix

    M = np.trace(np.dot(coil,sigma)) # Detect
    M_list.append(M) # Append to FID array


M = np.array(M_list)

figure('Excitation Profile')
title('Excitation Profile')
plot(omega_array, np.real(M), label = 'real')
plot(omega_array, np.imag(M), label = 'imag')
plot(omega_array, np.abs(M), label = 'abs')
legend()
xlabel('Frequency')

#figure('spec')
#title('Spectrum')
#plot(f,np.real(spec))
#xlabel('Frequency (Hz)')

Px = expm(1j*(np.pi/2)*sigma_x)
Py = expm(1j*(np.pi/2)*sigma_y)

#sigma_init = np.r_[
#        [
#            [1,0], 
#            [0,0]]
#        ]

sigma_init = sigma_z
#testx = Px @ sigma_init @ Px.conj()
#testy = Py @ sigma_init @ Py.conj()
testx = Px @ sigma_init @ Px.T
testy = Py @ sigma_init @ Py.T


print(testx)
print(testy)

show()