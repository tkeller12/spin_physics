import numpy as np
import matplotlib.pyplot as plt
from qutip import sigmax, sigmay, sigmaz, mesolve, Qobj, destroy

# Define system parameters
delta = np.linspace(-100e6, 100e6, 128)  # Detuning as an array of 10 points
omega_max = 2 * np.pi  # Maximum Rabi frequency in Hz
N = 1024  # Number of time steps
T1 = 1e-3  # Relaxation time in seconds
T2 = 5e-6  # Dephasing time in seconds
T = 10e-6  # Total simulation duration in seconds

tp = 100e-9
B1 = np.pi / tp

#def shaped_pulse(t, args):
#    """Gaussian envelope for the pulse."""
#    sigma = T / 20  # Standard deviation of the Gaussian pulse
#    return omega_max * np.exp(-0.5 * ((t - T/2) / sigma) ** 2)
def shaped_pulse(t, args):
#    """Gaussian envelope for the pulse."""
#    sigma = T / 20  # Standard deviation of the Gaussian pulse
    return B1 * int(t < tp)

# Define initial state
psi0 = Qobj([[1], [0]])  # Spin-up state

# Define collapse operators for relaxation and dephasing
c_ops = []
if T1 > 0:
    c_ops.append(np.sqrt(1/T1) * destroy(2))  # Relaxation
if T2 > 0:
    c_ops.append(np.sqrt(1/T2) * sigmaz())  # Dephasing

# Define time evolution
times = np.linspace(0, T, N)
results = []
for d in delta:
    H0 = d * sigmaz() / 2  # Static Hamiltonian
    H1 = sigmax() / 2  # Drive term
    H = [H0, [H1, shaped_pulse]]
    result = mesolve(H, psi0, times, c_ops, [sigmax(), sigmay(), sigmaz()])
    results.append(result.expect)

results = np.array(results)
Mx = results[:,0,:]
My = results[:,1,:]
Mz = results[:,2,:]

plt.figure()
plt.imshow(Mz, aspect = 'auto')
plt.colorbar()
plt.figure()
plt.imshow(My, aspect = 'auto')
plt.colorbar()

print(np.shape(Mz))
plt.figure()
#for ix, t_value in enumerate(times):
plt.plot(delta, Mx[:,-1], label = 'Mx')
plt.plot(delta, My[:,-1], label = 'My')
plt.plot(delta, Mz[:,-1], label = 'Mz')



## Plot results
#plt.figure(figsize=(10, 6))
#for i, d in enumerate(delta):
#    plt.plot(times, results[i][2], label=f'? = {d:.2f}')
#plt.xlabel('Time (s)')
#plt.ylabel('Expectation value of Z')
#plt.legend()
#plt.title('Spin-1/2 Dynamics under Shaped Pulse with Relaxation')
plt.show()
