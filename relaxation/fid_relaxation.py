import numpy as np
import matplotlib.pyplot as plt
from qutip import *

# ---------------------------
# Simulation parameters
# ---------------------------
T1 = 2.0   # longitudinal (T1) relaxation time in seconds
T2 = 1.0   # transverse (T2) relaxation time in seconds

# For a two-level system, 1/T2 = 1/(2T1) + 1/Tphi.
# Solve for Tphi:
Tphi = 1.0 / (1.0/T2 - 1.0/(2.0*T1))
print("Calculated pure dephasing time Tphi =", Tphi)

omega0 = 2 * np.pi   # Larmor frequency (rad/s)

# ---------------------------
# System definition
# ---------------------------
# Hamiltonian for a spin-1/2 in a static magnetic field along z:
H = 0.5 * omega0 * sigmaz()

print('H:', H)

# Initial state: after a 90° pulse, assume magnetization along x:
# ρ0 = 1/2*(I + σx)
rho0 = 0.5 * (qeye(2) + sigmax())

# ---------------------------
# Define relaxation (collapse) operators
# ---------------------------
c_ops = []
# T1 (energy relaxation): use sigma_minus with rate sqrt(1/T1)
c_ops.append(np.sqrt(1.0/T1) * sigmam())

# Pure dephasing: The Lindblad operator for pure dephasing is typically chosen as
# L_deph = sqrt(gamma_phi) * sigmaz, and it damps off-diagonals at a rate 2*gamma_phi.
# We require 2*gamma_phi = 1/Tphi, so gamma_phi = 1/(2*Tphi).
c_ops.append(np.sqrt(1.0/(2*Tphi)) * sigmaz())

# ---------------------------
# Time evolution parameters
# ---------------------------
tlist = np.linspace(0, 5, 200)  # simulate from 0 to 5 seconds

# We will monitor the expectation values of sigma_x, sigma_y, and sigma_z.
e_ops = [sigmax(), sigmay(), sigmaz()]

# ---------------------------
# Solve the master equation
# ---------------------------
result = mesolve(H, rho0, tlist, c_ops, e_ops)

# ---------------------------
# Plot the FID (transverse magnetization)
# ---------------------------
plt.figure(figsize=(8, 4))
plt.plot(tlist, result.expect[0], label=r'$\langle\sigma_x\rangle$')
plt.plot(tlist, result.expect[1], label=r'$\langle\sigma_y\rangle$')
plt.xlabel('Time (s)')
plt.ylabel('Magnetization')
plt.title('Simulated FID with Relaxation')
plt.legend()
plt.show()
