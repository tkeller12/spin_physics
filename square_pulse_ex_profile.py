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

print('Density Operators:')
print(sigma_x)
print(sigma_y)
print(sigma_z)
print('')
print('Pass Commutator Test:')
print(np.allclose(sigma_z * 1j, np.dot(sigma_x,sigma_y) - np.dot(sigma_y,sigma_x)))
print('')

omega = 500e6 # Bandwidth of Simulation (offset from carrier, Hz)
tp = 20e-9 # Pulse Length, s
pts = 1024 # Points in FID

pulse_angle = np.pi # pi/2 for 90-pulse, pi for 180-pulse

B1 = pulse_angle / tp

omega_array = np.r_[-omega/2:omega/2:1j*pts]

coil = sigma_x + 1j*sigma_y # Detection Operator (NMR Coil)

M_list = []
Mz_list = []

for ix in range(pts):
    sigma = sigma_z # Initial Density Matrix

    # re-calculate spin hamiltonian for offset
    H = 2*np.pi * tp * omega_array[ix] * sigma_z + B1 * tp * sigma_x # Calculate Hamiltonian (only Zeeman)
    P = expm(1j*H) # Define Propagator
    sigma = np.dot(np.dot(P,sigma),P.T.conj()) # Propagate Density Matrix

    M = np.trace(np.dot(coil,sigma)) # Detect
    Mz = np.trace(np.dot(sigma_z,sigma))
    M_list.append(M) # Append to FID array
    Mz_list.append(Mz)


M_array = np.array(M_list)
Mz_array = np.array(Mz_list)

figure('Excitation Profile')
title('Excitation Profile')
plot(omega_array/1e6, np.real(M_array), label = 'Mx')
plot(omega_array/1e6, np.imag(M_array), label = 'My')
plot(omega_array/1e6, np.real(Mz_array), label = 'Mz')
legend()
xlabel('Frequency Offset (MHz)')
show()
