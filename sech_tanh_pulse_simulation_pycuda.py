#import numpy as np
import cupy as cp
#from scipy.linalg import expm
from cupyx.scipy.linalg import expm
from matplotlib.pylab import *
import time

start_time = time.time()

sigma_x = 0.5*cp.r_[
        [
            [0, 1], 
            [1, 0] ]
        ] + 0j

sigma_y = 0.5*cp.r_[
        [
            [0,-1j], 
            [1j, 0] ]
        ] + 0j

sigma_z = 0.5*cp.r_[
        [
            [1, 0], 
            [0, -1] ]
        ] + 0j

print(sigma_x)
print(sigma_y)
print(sigma_z)


print('commutator test')
print(cp.dot(sigma_x,sigma_y) - cp.dot(sigma_y,sigma_x))
print(sigma_z * 1j)
print(cp.allclose(sigma_z * 1j, cp.dot(sigma_x,sigma_y) - cp.dot(sigma_y,sigma_x)))

pts = 128 # Points in frequency domain sweep
#pts = 4096 # Points in frequency domain sweep
omega_bw = 500e6 # Bandwidth of Simulation, Hz

#omega_array = cp.r_[-omega_bw/2.0:omega_bw/2.0:1j*pts]
omega_array = cp.linspace(-omega_bw/2.0,omega_bw/2.0,pts)

coil = sigma_x + 1j*sigma_y # Detection Operator (NMR Coil)

sech = lambda x: 1./cp.cosh(x)

dt = 1.e-9 # pulse resolution, s
tp = 128e-9 # Pulse length, s
amp = .1 # Pulse amplitude (AWG output level from 0 to 1)
BW = 100e6 # Pulse Bandwidth (FWHM), Hz
beta = 10. # Pulse truncation parameter

B1_max = 50.e6 # Maximum MW Field Strength, Hz

beta_tp = float(beta)/tp
mu = cp.pi*BW/beta_tp

#t = cp.r_[0.:tp:dt]
t = cp.linspace(0.,tp, int(tp/dt))

pulse = amp*(sech(beta_tp*(t-0.5*tp)))**(1.+1.j*mu)

#figure()
#plot(t*1e9, cp.real(pulse), label = 'real')
#plot(t*1e9, cp.imag(pulse), label = 'imag')
#legend()
#xlabel('Time (ns)')
#ylabel('B1')
#tight_layout()
#grid(linestyle = ':')

pulse *= B1_max*2.*cp.pi # Amplify Pulse

M_list = []
Mz_list = []

Mz_array = cp.zeros((len(t),len(omega_array)))
Mx_array = cp.zeros((len(t),len(omega_array)))
My_array = cp.zeros((len(t),len(omega_array)))

for omega_ix,omega in enumerate(omega_array):
    print('Offset: %i of %i'%((omega_ix+1),len(omega_array)))

    sigma = sigma_z # Initial Density Matrix

    # re-calculate spin hamiltonian for offset
    for time_ix,time_value in enumerate(t):
        B1 = pulse[time_ix]
        H = 2*cp.pi * omega * sigma_z + cp.real(B1) * sigma_x + cp.imag(B1) * sigma_y # Calculate Hamiltonian (only Zeeman)
        P = expm(1j*H*dt) # Define Propagator
        sigma = cp.dot(cp.dot(P,sigma),P.T.conj()) # Propagate Density Matrix

        Mx_value = cp.real(cp.trace(cp.dot(sigma_x,sigma)))
        My_value = cp.real(cp.trace(cp.dot(sigma_y,sigma)))
        Mx_array[time_ix, omega_ix] = Mx_value
        My_array[time_ix, omega_ix] = My_value
        Mz_value = cp.real(cp.trace(cp.dot(sigma_z,sigma)))
        Mz_array[time_ix, omega_ix] = Mz_value

    M = cp.trace(cp.dot(coil,sigma)) # Detect
    Mz_value = cp.trace(cp.dot(sigma_z,sigma))
    M_list.append(M) # Append to FID array
    Mz_list.append(Mz_value)

M = cp.array(M_list)
Mz = cp.array(Mz_list)
#
#figure('Excitation_Profile_Adiabatic_B1%0.0fMHz_%0.0fns_%0.0fMHz'%(B1_max/1e6,tp*1e9,BW/1e6))
#title('Excitation Profile, Adiabatic Pulse\n$B_{1}$=%0.0f MHz, tp=%0.0f ns, BW=%0.0f MHz, $\\beta$=%0.0f'%(B1_max/1e6,tp*1e9,BW/1e6, beta))
#plot(omega_array/1e6, cp.real(M), label = 'Mx')
#plot(omega_array/1e6, cp.imag(M), label = 'My')
#plot(omega_array/1e6, cp.real(Mz), label = 'Mz')
#legend()
#xlabel('Frequency (MHz)')
#
#figure()
#title('Mz')
#imshow(cp.real(Mz_array), aspect = 'auto', extent = (-omega_bw/2./1e6,omega_bw/2./1e6,tp*1e9,0.))
#xlabel('Frequency (MHz)')
#ylabel('Time (ns)')
#colorbar()
#
#ix = int(pts/2)
#freq = omega_array[ix]
#ax = figure().add_subplot(projection='3d')
#title('Offset %0.0f MHz'%(freq/1e6))
#ax.plot(cp.real(Mx_array)[:,ix], cp.real(My_array)[:,ix], cp.real(Mz_array)[:,ix])
#ax.set_xlim(-.5,.5)
#ax.set_ylim(-.5,.5)
#ax.set_zlim(-.5,.5)
#ax.set_xlabel('Mx')
#ax.set_ylabel('My')
#ax.set_zlabel('Mz')
#tight_layout()
stop_time = time.time()
print('Time Elapsed: %0.03f'%(stop_time - start_time))
show()
