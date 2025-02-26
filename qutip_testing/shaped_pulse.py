import numpy as np
import matplotlib.pyplot as plt
from qutip import sigmax, sigmay, sigmaz, mesolve, Qobj
import sys
import time
start_time = time.time()

# Define system parameters
delta_pts = 128
delta = np.linspace(-20, 20, delta_pts)  # Detuning as an array of 10 points
omega_max = 2 * np.pi * 5  # Maximum Rabi frequency in Hz
T = 1.0  # Total pulse duration in seconds
N = 128  # Number of time steps

#def shaped_pulse(t, args):
def shaped_pulse(t):
    """Gaussian envelope for the pulse."""
    sigma = T / 5  # Standard deviation of the Gaussian pulse
    return omega_max * np.exp(-0.5 * ((t - T/2) / sigma) ** 2)

# Define initial state
psi0 = Qobj([[1], [0]])  # Spin-up state

# Define time evolution
times = np.linspace(0, T, N)
results = []
last_updated = time.time()
for ix, d in enumerate(delta):
    if (time.time() - last_updated) > 0.2:
        print('%i of %i'%(ix, delta_pts))
        last_updated = time.time()
    H0 = d * sigmaz() / 2  # Static Hamiltonian
    H1 = sigmax() / 2  # Drive term
    H = [H0, [H1, shaped_pulse]]
    result = mesolve(H, psi0, times, [], [sigmax(), sigmay(), sigmaz()])
    results.append(result.expect)

# Plot results
#plt.figure(figsize=(10, 6))
#print(results)

results = np.array(results)
Mz = results[:,2,:]

plt.figure()
#plt.imshow(results[:][2][:], aspect = 'auto')
plt.imshow(Mz, aspect = 'auto')
#for i, d in enumerate(delta):
#    plt.plot(times, results[i][2], label=f'? = {d:.2f}')
#plt.xlabel('Time (s)')
#plt.ylabel('Expectation value of Z')
#plt.legend()
#plt.title('Spin-1/2 Dynamics under Shaped Pulse for Different Detunings')
stop_time = time.time()
print('Time Elapsed: %0.3f'%(stop_time - start_time))
plt.show()
