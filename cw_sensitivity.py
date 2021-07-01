import numpy as np

N = 1e19# Number of spins per unit volume
k_B = 1.38064852e-23 # Boltzmann's Constant, m^2 kg s^-2 K^-1
gamma = 1.76085964411e11 # gyromagnetic ratio, rad s^-1 T^-1
planck_const = 6.62607e-34 # Planck's constant
T = 300 # Absolute Temperature, K
frequency = 200e6# frequency
linewidth = 10e6# Lorentzian linewidth half height
filling_factor = 0.25 # filling factor
Q = 30 # Loaded Q-factor
Z0 = 50 # Characteristic impedance, ohms
P = 0.001 # Power, W

V = (N * (gamma**2.) * (planck_const / (2. * np.pi))**2.) / (4. * k_B * T) * (frequency / linewidth) * filling_factor * Q * np.sqrt(Z0*P)

print('Signal Voltage: %0.03f nV'%(V*1e9))
