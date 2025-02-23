import numpy as np
import matplotlib.pyplot as plt
from qutip import *

# ---------------------------
# Simulation parameters
# ---------------------------
T1 = 2.0       # T1 relaxation time (seconds)
T2 = 1.0       # T2 relaxation time (seconds)
# For a two-level system: 1/T2 = 1/(2T1) + 1/Tphi  =>  Tphi = 1/(1/T2 - 1/(2T1))
Tphi = 1.0 / (1.0/T2 - 1.0/(2.0*T1))
print("Pure dephasing time Tphi =", Tphi)

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
tlist = np.linspace(0, 5, 200)  # time from 0 to 5 seconds

# We will compute the expectation values of the transverse magnetizations for each spin.
# For spin 1, use: tensor(sx, I2); for spin 2: tensor(I2, sx).
e_ops = [tensor(sx, I2), tensor(I2, sx), tensor(sy, I2), tensor(I2, sy)]

# ---------------------------
# Solve the master equation using mesolve
# ---------------------------
result = mesolve(H, rho0, tlist, c_ops, e_ops)

# ---------------------------
# Plot the FID signal: transverse magnetization from both spins
# ---------------------------
plt.figure(figsize=(8, 4))
#plt.plot(tlist, result.expect[0], label=r'$\langle\sigma_x^{(1)}\rangle$')
#plt.plot(tlist, result.expect[1], label=r'$\langle\sigma_x^{(2)}\rangle$')
plt.plot(tlist, result.expect[0]+result.expect[1], label=r'$\langle\sigma_x^{(1)}\rangle+\langle\sigma_x^{(2)}\rangle$')
plt.plot(tlist, result.expect[2]+result.expect[3], label=r'$\langle\sigma_y^{(1)}\rangle+\langle\sigma_y^{(2)}\rangle$')
#plt.plot(tlist, result.expect[1], label=r'$\langle\sigma_x^{(2)}\rangle$')
plt.xlabel('Time (s)')
plt.ylabel('Transverse Magnetization')
plt.title('Simulated FID with J-Coupling and Relaxation (Two Spins)')
plt.legend()
plt.show()
