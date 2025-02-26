import concurrent.futures

import numpy as np
from scipy.linalg import expm
#from scipy.linalg import expm_cond as expm
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


pts = 16 # Points in frequency domain sweep
#pts = 4096 # Points in frequency domain sweep
omega_bw = 500e6 # Bandwidth of Simulation, Hz

omega_array = np.r_[-omega_bw/2:omega_bw/2:1j*pts]

coil = sigma_x + 1j*sigma_y # Detection Operator (NMR Coil)

sech = lambda x: 1./np.cosh(x)

dt = 1.e-9 # pulse resolution, s
tp = 128e-9 # Pulse length, s
amp = .1 # Pulse amplitude (AWG output level from 0 to 1)
BW = 100e6 # Pulse Bandwidth (FWHM), Hz
beta = 10. # Pulse truncation parameter

B1_max = 50.e6 # Maximum MW Field Strength, Hz

beta_tp = float(beta)/tp
mu = np.pi*BW/beta_tp

t = np.r_[0.:tp:dt]

pulse = amp*(sech(beta_tp*(t-0.5*tp)))**(1.+1.j*mu)

pulse *= B1_max*2.*np.pi # Amplify Pulse


def generate_array(omega_ix):
    sigma = sigma_z # Initial Density Matrix
    Mx = np.zeros_like(t)
    My = np.zeros_like(t)
    Mz = np.zeros_like(t)
    omega = omega_array[omega_ix]

    B1 = pulse[0]
    H = 2*np.pi * omega * sigma_z + np.real(B1) * sigma_x + np.imag(B1) * sigma_y # Calculate Hamiltonian (only Zeeman)
    P = expm(1j*H*dt) # Define Propagator
    # re-calculate spin hamiltonian for offset
    for time_ix,time in enumerate(t):
        B1 = pulse[time_ix]
        H = 2*np.pi * omega * sigma_z + np.real(B1) * sigma_x + np.imag(B1) * sigma_y # Calculate Hamiltonian (only Zeeman)
#        P = expm(1j*H*dt) # Define Propagator
        sigma = np.dot(np.dot(P,sigma),P.T.conj()) # Propagate Density Matrix

#        Mx_value = np.real(np.trace(np.dot(sigma_x,sigma)))
#        My_value = np.real(np.trace(np.dot(sigma_y,sigma)))
#        Mz_value = np.real(np.trace(np.dot(sigma_z,sigma)))
#        Mx[time_ix] = Mx_value
#        My[time_ix] = My_value
#        Mz[time_ix] = Mz_value

#    M = np.trace(np.dot(coil,sigma)) # Detect
#    Mz_value = np.trace(np.dot(sigma_z,sigma))

#    return omega_ix, Mx, My, Mz

#for omega_ix, omega in enumerate(omega_array):
#    print('Offset: %i of %i'%((omega_ix+1),len(omega_array)))
#    Mx, My, Mz = generate_array(omega_ix, omega)
#    Mx_array[:, omega_ix] = Mx
#    My_array[:, omega_ix] = My
#    Mz_array[:, omega_ix] = Mz


#figure('Excitation_Profile_Adiabatic_B1%0.0fMHz_%0.0fns_%0.0fMHz'%(B1_max/1e6,tp*1e9,BW/1e6))
#title('Excitation Profile, Adiabatic Pulse\n$B_{1}$=%0.0f MHz, tp=%0.0f ns, BW=%0.0f MHz, $\\beta$=%0.0f'%(B1_max/1e6,tp*1e9,BW/1e6, beta))
#plot(omega_array/1e6, np.real(M), label = 'Mx')
#plot(omega_array/1e6, np.imag(M), label = 'My')
#plot(omega_array/1e6, np.real(Mz), label = 'Mz')
#legend()
#xlabel('Frequency (MHz)')

def main():

#    array_size = pts  # Size of each generated array
#    num_tasks = 128    # Number of parallel tasks (and thus arrays)
#    seeds = [100 + i for i in range(num_tasks)]  # Unique seeds for each task
    
    # Initialize an array to accumulate the sum.
#    sum_array = np.zeros(array_size, dtype = np.complex64)
    num_finished = 0
    num_tasks = len(omega_array)
    
    # Use ProcessPoolExecutor for parallel processing.
    with concurrent.futures.ProcessPoolExecutor() as executor:
        # Process each generated array as soon as it's available.
        for each in executor.map(generate_array, range(num_tasks)):
            num_finished += 1
#            sum_array += arr
            print('%0.01f%%'%(100*num_finished / num_tasks))
#            Mx_array[:, omega_ix] = Mx
#            My_array[:, omega_ix] = My
#            Mz_array[:, omega_ix] = Mz
    
#    print("Sum of arrays:")
#    print(sum_array)

if __name__ == "__main__":
    Mz_array = np.zeros((len(t),len(omega_array)))
    Mx_array = np.zeros((len(t),len(omega_array)))
    My_array = np.zeros((len(t),len(omega_array)))
    start_time = time.time()
    main()
    stop_time = time.time()
    print('Time Elapsed: %0.03f s'%(stop_time - start_time))

#    figure()
#    title('Mz')
#    imshow(np.real(Mz_array), aspect = 'auto', extent = (-omega_bw/2./1e6,omega_bw/2./1e6,tp*1e9,0.))
#    xlabel('Frequency (MHz)')
#    ylabel('Time (ns)')
#    colorbar()
#    show()
#
#ix = int(pts/2)
#freq = omega_array[ix]
#ax = figure().add_subplot(projection='3d')
#title('Offset %0.0f MHz'%(freq/1e6))
#ax.plot(np.real(Mx_array)[:,ix], np.real(My_array)[:,ix], np.real(Mz_array)[:,ix])
#ax.set_xlim(-.5,.5)
#ax.set_ylim(-.5,.5)
#ax.set_zlim(-.5,.5)
#ax.set_xlabel('Mx')
#ax.set_ylabel('My')
#ax.set_zlabel('Mz')
#tight_layout()
