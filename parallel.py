import concurrent.futures

import numpy as np
import matplotlib.pyplot as plt
from qutip import *

print('Initializing parameters...')
# ---------------------------
# Simulation parameters
# ---------------------------
T1 = 2.0       # T1 relaxation time (seconds)
T2 = 1.0       # T2 relaxation time (seconds)
# For a two-level system: 1/T2 = 1/(2T1) + 1/Tphi  =>  Tphi = 1/(1/T2 - 1/(2T1))
Tphi = 1.0 / (1.0/T2 - 1.0/(2.0*T1))
#print("Pure dephasing time Tphi =", Tphi)

omega0 = 2 * np.pi   # Larmor frequency (rad/s)
J = 1.1              # J-coupling in Hz (for simplicity)

# ---------------------------
# Define operators for a single spin-1/2
# ---------------------------
sx = sigmax()
sy = sigmay()
sz = sigmaz()
sm = sigmam()
I2 = qeye(2)

# ---------------------------
# Hamiltonian for two coupled spins
# ---------------------------
# Zeeman terms for spin 1 and spin 2:
H1 = 0.5 * omega0 * tensor(sz, I2)
H2 = 0.5 * (omega0*1.1) * tensor(I2, sz)

# J-coupling: H_J = 2*pi*J*I_{1z} I_{2z} = (2*pi*J/4) * tensor(sz, sz)
HJ = (2 * np.pi * J / 4.0) * tensor(sz, sz)

# Total Hamiltonian:
H = H1 + H2 + HJ

# ---------------------------
# Initial state: Both spins rotated 90° into the x-axis.
# For a single spin: rho = 1/2*(I + sigma_x). For two spins, take the tensor product.
# ---------------------------
rho0_spin = 0.5 * (I2 + sx)
rho0 = tensor(rho0_spin, rho0_spin)

# ---------------------------
# Collapse operators for relaxation
# ---------------------------
c_ops = []
# T1 relaxation for spin 1 and spin 2:
c_ops.append(np.sqrt(1.0/T1) * tensor(sm, I2))
c_ops.append(np.sqrt(1.0/T1) * tensor(I2, sm))
# Pure dephasing for spin 1 and spin 2:
c_ops.append(np.sqrt(1.0/(2*Tphi)) * tensor(sz, I2))
c_ops.append(np.sqrt(1.0/(2*Tphi)) * tensor(I2, sz))

# ---------------------------
# Time evolution parameters
# ---------------------------
pts = 4096
tlist = np.linspace(0, 5, pts)  # time from 0 to 5 seconds

# We will compute the expectation values of the transverse magnetizations for each spin.
# For spin 1, use: tensor(sx, I2); for spin 2: tensor(I2, sx).
e_ops = [tensor(sx, I2), tensor(I2, sx), tensor(sy, I2), tensor(I2, sy)]

print('Done.')



def generate_array(size = 100, seed = 0):
    """
    Generate a NumPy array of random numbers with the given size.

    Args:
        size (int): Size of the 1D NumPy array.
        seed (int): Seed for the random number generator for reproducibility.

    Returns:
        np.ndarray: A 1D array of random numbers.
    """
    result = mesolve(H, rho0, tlist, c_ops, e_ops = e_ops)
    out = result.expect[0]+result.expect[1] + 1j * (result.expect[2]+result.expect[3])
    return np.ones_like(out)

def main():

    array_size = pts  # Size of each generated array
    num_tasks = 128    # Number of parallel tasks (and thus arrays)
    seeds = [100 + i for i in range(num_tasks)]  # Unique seeds for each task
    
    # Initialize an array to accumulate the sum.
    sum_array = np.zeros(array_size, dtype = np.complex64)
    num_finished = 0
    
    # Use ProcessPoolExecutor for parallel processing.
    with concurrent.futures.ProcessPoolExecutor() as executor:
        # Process each generated array as soon as it's available.
        for arr in executor.map(generate_array, [array_size] * num_tasks, seeds):
            num_finished += 1
            sum_array += arr
            print('%0.01f%%'%(100*num_finished / num_tasks))
    
    print("Sum of arrays:")
    print(sum_array)

if __name__ == "__main__":
    main()
