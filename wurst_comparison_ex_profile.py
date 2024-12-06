import numpy as np
from scipy.linalg import expm
from matplotlib.pylab import *

def wurst(tp, N, resolution):
    '''Real value WURST envelope pulse shape

    .. math::
        1 - \\text{abs} \left( \cos \left( \\frac{\pi}{t_p} (t - \\frac{t_p}{2}) + \\frac{\pi}{2} \\right) \\right) ^N

    Args:
        tp (float): Pulse length
        N (float): exponential 

    Returns:
        tuple: tuple containing:

        t (*numpy.ndarray*): Time axis

        pulse (*numpy.ndarray*): Pulse shape
    '''
    t = np.r_[0.:tp:resolution]
    pulse = (1. - np.abs(np.cos(np.pi*(t-tp/2.)/tp + np.pi/2.))**N) + 0j

    return t, pulse

def chirp(tp, BW, resolution):
    '''Complex chirp pulse

    .. math::
        e^{i 2 \pi (k/2) (t - t_p/2)^2}

    Args:
        tp (float): Pulse length
        BW (float): Bandwidth of pulse

    Returns:
        tuple: tuple containing:

            t (*numpy.ndarray*): Time axis

            pulse (*numpy.ndarray*): Pulse shape
    '''
    k = BW/tp
    t = np.r_[0.:tp:resolution]
    pulse = np.exp(1.j*2.*np.pi*((k/2.)*((t-tp/2.)**2.)))
    return t, pulse


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

print('Sx: ', sigma_x)
print('Sy: ', sigma_y)
print('Sz: ', sigma_z)


print('Test Commutator:')
print(np.dot(sigma_x,sigma_y) - np.dot(sigma_y,sigma_x))
print(sigma_z * 1j)
print(np.allclose(sigma_z * 1j, np.dot(sigma_x,sigma_y) - np.dot(sigma_y,sigma_x)))

omega_bw = 500e6 # frequency offset from carrier, Hz

tp_180 = 10e-9# 180-pulse length, s
pts = 128 # Points in FID

pulse_B1 = np.pi / tp_180
print('B1 power: %0.03e Angular Frequency'%pulse_B1)

omega_array = np.r_[-omega_bw/2:omega_bw/2:1j*pts]

coil = sigma_x + 1j*sigma_y # Detection Operator (NMR Coil)

for power in [5, 10]:
    print('Power: ', power)

    tp = 128e-9
    BW = 220e6
    dt = 0.1e-9
    amp = 1.
    t, shape = wurst(tp, power, dt)
    t, chirp_shape = chirp(tp, BW, dt)
    pulse = pulse_B1 * amp * shape * chirp_shape

    figure('Pulse Shape')
    plot(t*1e9,np.real(pulse)/pulse_B1, label = 'WURST-%i'%power)
    xlabel('Time (ns)')
    ylabel('Normalized Pulse Amplitude')
    legend()
    xlim(np.min(t*1e9), np.max(t*1e9))
    tight_layout()

    M_list = []
    Mz_list = []

    Mz_array = np.zeros((len(t),len(omega_array)))

    for omega_ix,omega in enumerate(omega_array):
        sigma = sigma_z # Initial Density Matrix

        for time_ix,time in enumerate(t):
            B1 = pulse[time_ix]
            # re-calculate spin hamiltonian for offset
            H = 2*np.pi * omega * sigma_z + np.real(B1) * sigma_x + np.imag(B1) * sigma_y # Calculate Hamiltonian
            P = expm(1j*H*dt) # Define Propagator
            sigma = np.dot(np.dot(P,sigma),P.T.conj()) # Propagate Density Matrix

            Mz_value = 2*np.real(np.trace(np.dot(sigma_z,sigma)))
            Mz_array[time_ix, omega_ix] = Mz_value

        M = 2*np.trace(np.dot(coil,sigma)) # Detect Sx & Sy
    #    Mz_value = 2*np.trace(np.dot(sigma_z,sigma)) # Detect Sz
        M_list.append(M)
        Mz_list.append(Mz_value)

    M = np.array(M_list)
    Mz = np.array(Mz_list)

    #figure('Excitation_Profile_WURST-%i_%0.0fns_%0.0fMHz'%(power,tp*1e9,BW/1e6))
    #title('Excitation Profile\nWURST-%i, %0.0f ns, %0.0f MHz'%(power,tp*1e9,BW/1e6))
    figure('Excitation_Profile')
    title('Excitation Profile')
    #plot(omega_array/1e6, np.real(M), label = 'Mx')
    #plot(omega_array/1e6, np.imag(M), label = 'My')
    plot(omega_array/1e6, np.real(Mz), label = 'WURST-%i'%power)
    xlim(np.min(omega_array/1e6), np.max(omega_array/1e6))
    legend()
    xlabel('Frequency (MHz)')
    ylabel('Mz')
    tight_layout()

    #figure()
    #imshow(np.real(Mz_array), aspect = 'auto')
    #colorbar()
show()
