import numpy as np
from scipy.linalg import expm
from matplotlib.pylab import *

import pyDEER as deer

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

omega_bw = 500e6 # frequency offset from carrier, Hz
#tp = 0.0001 # Pulse Length, s
B1 = 125e6
#B1 = 0e6
pts = 128 # Points in FID

#tp = np.pi/2 / B1

omega_array = np.r_[-omega_bw/2:omega_bw/2:1j*pts]

coil = sigma_x + 1j*sigma_y # Detection Operator (NMR Coil)


tp = 200e-9
BW = 100e6
dt = 1e-9
amp = 1.
t,shape = deer.wurst(tp, 100, resolution = dt)
#t,shape = deer.adiabatic(tp, BW, 3, resolution = dt)
#t, shape = deer.sinc(tp, 10, resolution = dt)
t, chirp = deer.chirp(tp, BW, resolution = dt)

pulse = amp * shape * chirp
#pulse = amp * shape

pulse *= B1
figure('Pulse Shape')
plot(t,np.real(pulse))
plot(t,np.imag(pulse))
#H = 2*np.pi * omega * sigma_z # Calculate Hamiltonian (only Zeeman)

#P = expm(1j*H*dt) # Define Propagator

#sigma = sigma_z # Initial Density Matrix (After 90-pulse)

M_list = []
Mz_list = []

Mz_array = np.zeros((len(t),len(omega_array)))

for omega_ix,omega in enumerate(omega_array):
#    sigma = sigma_z # Initial Density Matrix

    sigma = np.r_[[[1,0], [0,0]]] + np.eye(2)
    # re-calculate spin hamiltonian for offset
#    H = 2*np.pi * omega_array[ix] * sigma_z + np.pi/2 * sigma_y # Calculate Hamiltonian (only Zeeman)

    for time_ix,time in enumerate(t):
        B1 = pulse[time_ix]
        H = 2*np.pi * omega * sigma_z + np.real(B1) * sigma_x + np.imag(B1) * sigma_y # Calculate Hamiltonian (only Zeeman)
        P = expm(1j*H*dt) # Define Propagator
        sigma = np.dot(np.dot(P,sigma),P.T.conj()) # Propagate Density Matrix

        Mz_value = np.real(np.trace(np.dot(sigma_z,sigma)))
        Mz_array[time_ix, omega_ix] = Mz_value
#    sigma = P @ sigma @ P.T.conj()

    M = np.trace(np.dot(coil,sigma)) # Detect
    Mz_value = np.trace(np.dot(sigma_z,sigma))
    M_list.append(M) # Append to FID array
    Mz_list.append(Mz_value)


M = np.array(M_list)
Mz = np.array(Mz_list)

figure('Excitation Profile')
title('Excitation Profile')
plot(omega_array, np.real(M), label = 'Mx')
plot(omega_array, np.imag(M), label = 'My')
plot(omega_array, np.real(Mz), label = 'Mz')
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
#sigma_init = np.r_[[[1,0], [0,0]]]
#testx = Px @ sigma_init @ Px.conj()
#testy = Py @ sigma_init @ Py.conj()
testx = Px @ sigma_init @ Px.T # x-pulse does not work
testy = Py @ sigma_init @ Py.T


print(testx)
print(testy)

figure()
imshow(np.real(Mz_array), aspect = 'auto')
colorbar()


show()
