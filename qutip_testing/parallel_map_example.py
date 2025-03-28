import numpy as np
import qutip as qt
from qutip import parallel_map

# Define system parameters
N = 10
omega = 1.0
timesteps = np.linspace(0, 10, 100)

# Define Hamiltonian and initial state globally
H = omega * qt.sigmaz()
psi0 = qt.basis(2, 0)

def evolve_and_measure(t):
    """ Evolves the state psi0 under Hamiltonian H for time t and computes expectation value. """
    result = qt.sesolve(H, psi0, [t], [qt.sigmax()])
    return result.expect[0][0]

if __name__ == "__main__":  # Ensures safe multiprocessing
    results = parallel_map(evolve_and_measure, timesteps)
    print(results[:10])

