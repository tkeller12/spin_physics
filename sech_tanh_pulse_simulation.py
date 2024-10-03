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

pts = 128 # Points in frequency domain sweep
omega_bw = 500e6 # Bandwidth of Simulation, Hz

omega_array = np.r_[-omega_bw/2:omega_bw/2:1j*pts]

coil = sigma_x + 1j*sigma_y # Detection Operator (NMR Coil)


sech = lambda x: 1./np.cosh(x)

dt = 1.e-9 # pulse resolution, s
tp = 128e-9 # Pulse length, s
amp = 1. # Pulse amplitude (AWG output level from 0 to 1)
BW = 100e6 # Pulse Bandwidth (FWHM), Hz
beta = 10. # Pulse truncation parameter

B1_max = 50.e6 # Maximum MW Field Strength, Hz

beta_tp = float(beta)/tp
mu = np.pi*BW/beta_tp

t = np.r_[0.:tp:dt]

pulse = (sech(beta_tp*(t-0.5*tp)))**(1.+1.j*mu)

figure()
plot(t*1e9, np.real(pulse), label = 'real')
plot(t*1e9, np.imag(pulse), label = 'imag')
legend()
xlabel('Time (ns)')
ylabel('B1')
tight_layout()
grid(linestyle = ':')

pulse *= B1_max*2.*np.pi # Amplify Pulse

M_list = []
Mz_list = []

Mz_array = np.zeros((len(t),len(omega_array)))

for omega_ix,omega in enumerate(omega_array):
    print('Offset: %i of %i'%((omega_ix+1),len(omega_array)))

    sigma = sigma_z # Initial Density Matrix

    # re-calculate spin hamiltonian for offset
    for time_ix,time in enumerate(t):
        B1 = pulse[time_ix]
        H = 2*np.pi * omega * sigma_z + np.real(B1) * sigma_x + np.imag(B1) * sigma_y # Calculate Hamiltonian (only Zeeman)
        P = expm(1j*H*dt) # Define Propagator
        sigma = np.dot(np.dot(P,sigma),P.T.conj()) # Propagate Density Matrix

        Mz_value = np.real(np.trace(np.dot(sigma_z,sigma)))
        Mz_array[time_ix, omega_ix] = Mz_value

    M = np.trace(np.dot(coil,sigma)) # Detect
    Mz_value = np.trace(np.dot(sigma_z,sigma))
    M_list.append(M) # Append to FID array
    Mz_list.append(Mz_value)

M = np.array(M_list)
Mz = np.array(Mz_list)

figure('Excitation_Profile_Adiabatic_B1%0.0fMHz_%0.0fns_%0.0fMHz'%(B1_max/1e6,tp*1e9,BW/1e6))
title('Excitation Profile, Adiabatic Pulse\n$B_{1}$=%0.0f MHz, tp=%0.0f ns, BW=%0.0f MHz, $\\beta$=%0.0f'%(B1_max/1e6,tp*1e9,BW/1e6, beta))
plot(omega_array/1e6, np.real(M), label = 'Mx')
plot(omega_array/1e6, np.imag(M), label = 'My')
plot(omega_array/1e6, np.real(Mz), label = 'Mz')
legend()
xlabel('Frequency (MHz)')

figure()
title('Mz')
imshow(np.real(Mz_array), aspect = 'auto', extent = (-omega_bw/2./1e6,omega_bw/2./1e6,tp*1e9,0.))
xlabel('Frequency (MHz)')
ylabel('Time (ns)')
colorbar()

show()
