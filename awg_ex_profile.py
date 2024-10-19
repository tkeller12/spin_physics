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

omega = 200e6 # frequency offset from carrier, Hz
#tp = 0.0001 # Pulse Length, s
pts = 100 # Points in freq

#tp = np.pi/2 / B1

pulse = awg.gaussian_pulse(50e-9,5)
#pulse = awg.sinc(50e-9,9)
time = pulse[0]
shape = pulse[1]
B1 = 35e6

dt = time[1] - time[0]

figure('pulse shape')
plot(time,np.real(shape))
plot(time,np.imag(shape))
#show()

freq = np.r_[-omega/2:omega/2:1j*pts]

coil = sigma_x + 1j*sigma_y # Detection Operator (NMR Coil)

#H = 2*np.pi * omega * sigma_z # Calculate Hamiltonian (only Zeeman)

#P = expm(1j*H*dt) # Define Propagator

#sigma = sigma_z # Initial Density Matrix (After 90-pulse)

M = np.zeros((len(time),len(freq)))
Mz = np.zeros((len(time),len(freq)))

#M_list = []

for f_ix, f in enumerate(freq):
    print(f_ix)
    sigma = sigma_z # Initial Density Matrix

#    sigma = np.r_[[[1,0], [0,0]]]

    for t_ix, t in enumerate(time):
        # re-calculate spin hamiltonian for offset

        H = 2*np.pi * freq[f_ix] * sigma_z
        H += 2*np.pi * B1 * np.real(shape[t_ix]) * sigma_x
        H += 2*np.pi * B1 * np.imag(shape[t_ix]) * sigma_y

        P = expm(-1j*dt*H) # Define Propagator

        sigma =  P @ sigma @ P.T.conj()

        M[t_ix,f_ix] = np.trace(np.dot(coil,sigma)) # Detect
        Mz[t_ix,f_ix] = np.trace(np.dot(sigma_z,sigma)) # Detect
        

extent = np.r_[np.min(freq)*1e-6, np.max(freq)*1e-6, np.max(time)*1e9, np.min(time)*1e9]

figure('real M, Mx')
imshow(np.abs(M), aspect = 'auto', extent = extent)
xlabel('Frequency (MHz)')
ylabel('Time (ns)')
colorbar()

figure('imag M, My')
imshow(np.imag(M), aspect = 'auto', extent = extent)
xlabel('Frequency (MHz)')
ylabel('Time (ns)')
colorbar()

figure('Mz')
imshow(np.real(Mz), aspect = 'auto', extent = extent)
xlabel('Frequency (MHz)')
ylabel('Time (ns)')
colorbar()

M_end = M[-1,:].ravel()
Mz_end = Mz[-1,:].ravel()

figure()
plot(freq/1e6, np.real(M_end), label = 'Mx')
plot(freq/1e6, np.imag(M_end), label = 'My')
plot(freq/1e6, np.real(Mz_end), label = 'Mz')
xlabel('Frequency (Mz)')
ylabel('M')
legend()
show()

show()
