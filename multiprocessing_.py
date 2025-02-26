import multiprocessing as mp
import numpy as np
import math
import time

import numpy as np
from scipy.linalg import expm
from matplotlib.pylab import *
print('Starting Script...')

def compute_value(omega_ix):
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
    sigma = sigma_z # Initial Density Matrix
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

#    B1 = pulse[0]
    Mx = np.zeros_like(t)
    My = np.zeros_like(t)
    Mz = np.zeros_like(t)
    omega = omega_array[omega_ix]

    # re-calculate spin hamiltonian for offset
    for time_ix,time in enumerate(t):
        B1 = pulse[time_ix]
        H = 2*np.pi * omega * sigma_z + np.real(B1) * sigma_x + np.imag(B1) * sigma_y # Calculate Hamiltonian (only Zeeman)
        P = expm(1j*H*dt) # Define Propagator
        sigma = np.dot(np.dot(P,sigma),P.T.conj()) # Propagate Density Matrix

    return omega_ix

def main():
    # Define the size of the array to generate.
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

    num_elements = pts  # Adjust this value for a larger array/longer computation
    print('main')



    
    # Create a Pool with as many processes as there are CPU cores.
    with mp.Pool(processes=mp.cpu_count()) as pool:
        start_time = time.time()
        
        # Map the compute_value function to each index in parallel.
        results = pool.map(compute_value, range(num_elements))
        print(results)
        
        elapsed_time = time.time() - start_time
        print(f"Computation completed in {elapsed_time:.2f} seconds.")
    
    # Convert the list of results to a NumPy array.
    result_array = np.array(results)
    print("Result array:")
    print(result_array)

if __name__ == '__main__':
    main()
