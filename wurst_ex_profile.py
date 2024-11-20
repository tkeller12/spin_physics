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

tp_180 = 10e-9# 180-pulse length, s
pts = 128 # Points in FID

#tp = np.pi/2 / B1
B1 = np.pi / tp_180
print('B1 power: %0.03e Angular Frequency'%B1)

omega_array = np.r_[-omega_bw/2:omega_bw/2:1j*pts]

coil = sigma_x + 1j*sigma_y # Detection Operator (NMR Coil)


tp = 128e-9
BW = 100e6
dt = 0.5e-9
amp = 0.3
power = 5. # 
t,shape = deer.wurst(tp, power, resolution = dt)
#t,shape = deer.adiabatic(tp, BW, 3, resolution = dt)
#t, shape = deer.sinc(tp, 10, resolution = dt)
t, chirp = deer.chirp(tp, BW, resolution = dt)

pulse = amp * shape * chirp
#pulse = amp * chirp
#pulse = amp * shape

pulse *= B1
figure('Pulse Shape')
plot(t,np.real(pulse))
plot(t,np.imag(pulse))

M_list = []
Mz_list = []

Mz_array = np.zeros((len(t),len(omega_array)))

for omega_ix,omega in enumerate(omega_array):
    print('%i of %i'%(omega_ix,len(omega_array)))
    sigma = sigma_z # Initial Density Matrix

#    sigma = np.r_[[[1,0], [0,0]]] + np.eye(2)
    # re-calculate spin hamiltonian for offset

    for time_ix,time in enumerate(t):
        B1 = pulse[time_ix]
        H = 2*np.pi * omega * sigma_z + np.real(B1) * sigma_x + np.imag(B1) * sigma_y # Calculate Hamiltonian (only Zeeman)
        P = expm(1j*H*dt) # Define Propagator
        sigma = np.dot(np.dot(P,sigma),P.T.conj()) # Propagate Density Matrix

        Mz_value = np.real(np.trace(np.dot(sigma_z,sigma)))
        Mz_array[time_ix, omega_ix] = Mz_value

    M = 2*np.trace(np.dot(coil,sigma)) # Detect
    Mz_value = 2*np.trace(np.dot(sigma_z,sigma))
    M_list.append(M) # Append to FID array
    Mz_list.append(Mz_value)


M = np.array(M_list)
Mz = np.array(Mz_list)

figure('Excitation_Profile_WURST-%i_%0.0fns_%0.0fMHz'%(power,tp*1e9,BW/1e6))
title('Excitation Profile\nWURST-%i, %0.0f ns, %0.0f MHz'%(power,tp*1e9,BW/1e6))
plot(omega_array/1e6, np.real(M), label = 'Mx')
plot(omega_array/1e6, np.imag(M), label = 'My')
plot(omega_array/1e6, np.real(Mz), label = 'Mz')
legend()
xlabel('Frequency (MHz)')

figure()
imshow(np.real(Mz_array), aspect = 'auto')
colorbar()

show()
