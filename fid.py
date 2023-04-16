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
print(np.allclose(sigma_z * 1j, np.dot(sigma_x,sigma_y) - np.dot(sigma_y,sigma_x)))

omega = 2. # frequency offset from carrier, Hz
dt = 0.01 # Dwell Time, s
pts = 1024 # Points in FID
T2 = 2. # T2 used for apodization
t = np.r_[0:pts] * dt

coil = sigma_x + 1j*sigma_y # Detection Operator (NMR Coil)

H = 2*np.pi * omega * sigma_z # Calculate Hamiltonian (only Zeeman)

P = expm(1j*H*dt) # Define Propagator

sigma = sigma_x # Initial Density Matrix (After 90-pulse)

M_list = []

for ix in range(pts):
    M = np.trace(np.dot(coil,sigma)) # Detect
    M_list.append(M) # Append to FID array
    sigma = np.dot(np.dot(P,sigma),np.conjugate(P)) # Propagate Density Matrix


M = np.array(M_list)
fid = M*np.exp(-1.*t/T2) # Apply Apodization

f = np.r_[-0.5/dt:0.5/dt:1j*pts] # Calculate frequency 
print(f)

spec = np.fft.fftshift(np.fft.fft(fid)) # Fourier Transform

figure('fid')
title('FID')
plot(t, np.real(fid), label = 'real')
plot(t, np.imag(fid), label = 'imag')
legend()
xlabel('Time (s)')

figure('spec')
title('Spectrum')
plot(f,np.real(spec))
xlabel('Frequency (Hz)')
show()
