import numpy as np
import matplotlib.pyplot as plt
from qutip import sigmax, sigmay, sigmaz, mesolve, Qobj, destroy, parallel_map, serial_map

# Define system parameters
#delta = np.linspace(-100e6, 100e6, 4096)
#delta = np.linspace(-200e6, 200e6, 1025)
#delta = np.linspace(-200e6, 200e6, 1001)
delta = np.linspace(-250e6, 250e6, 1001)
#delta = np.linspace(-100e6, 100e6, 100)  # Detuning as an array of 10 points
omega_max = 2 * np.pi  # Maximum Rabi frequency in Hz
N = 1024  # Number of time steps
#N = 256  # Number of time steps
T1 = 10e-3  # Relaxation time in seconds
T2 = 2e-6  # Dephasing time in seconds

#T1 = 0e-3  # Relaxation time in seconds
#T2 = 0e-6  # Dephasing time in seconds
T = 10e-6  # Total simulation duration in seconds

tp = 100e-9
B1 = np.pi / tp

tau = 1000e-9
tau2 = 4000e-9

center = 0e6
width = 25e6
spectrum = np.exp(-0.5 * ((delta-center) / width)**2.0)

#def shaped_pulse(t, args):
#    """Gaussian envelope for the pulse."""
#    sigma = T / 20  # Standard deviation of the Gaussian pulse
#    return omega_max * np.exp(-0.5 * ((t - T/2) / sigma) ** 2)
#    """Gaussian envelope for the pulse."""
#    sigma = T / 20  # Standard deviation of the Gaussian pulse
#    return B1 * int(t < tp)
#    return B1 * np.round(t < tp)
#def shaped_pulse(t, args = None):
def shaped_pulse(t):
    p90_1 = B1/2 * np.array((t < tp), dtype = np.float64)
    p90_2 = B1/2 * np.array(np.abs(t - tau) < (tp/2), dtype = np.float64)
    p90_3 = B1/2 * np.array(np.abs(t - (tau+tau2)) < (tp/2), dtype = np.float64)
    return p90_1 + p90_2 + p90_3

# Define initial state
psi0 = Qobj([[1], [0]])  # Spin-up state

# Define collapse operators for relaxation and dephasing
c_ops = []
if T1 > 0:
#    c_ops.append(np.sqrt(1.0/T1) * destroy(2))  # Relaxation
    c_ops.append(np.sqrt(1.0/T1) * destroy(2))  # Relaxation
if T2 > 0:
#    c_ops.append(np.sqrt(1/T2) * sigmaz())  # Dephasing
    c_ops.append(np.sqrt(1.0/(2.0*T2)) * sigmaz())  # Dephasing

# Define time evolution
times = np.linspace(0, T, N)
pulse_test = shaped_pulse(times)
#print(times)
#print(pulse_test)
plt.figure('Pulse Shape')
plt.plot(times,pulse_test)
results = []

def sim(d):
    H0 = d * sigmaz() / 2.0  # Static Hamiltonian
    H1 = sigmax() / 2.0  # Drive term
    H = [H0, [H1, shaped_pulse]]
#    result = mesolve(H, psi0, times, c_ops, [sigmax(), sigmay(), sigmaz()])
    result = mesolve(H, psi0, times, c_ops, e_ops = [sigmax(), sigmay(), sigmaz()])
    return result.expect
#    return 0#result.expect
#    return result

#for d in delta:
#    result = sim(d)
#    results.append(result.expect)

if __name__ == "__main__":
    results = parallel_map(sim, delta)
    #results = serial_map(sim, delta)
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

#    print(np.shape(Mz))
    plt.figure('Mx, My, Mz')
    #for ix, t_value in enumerate(times):
#    ix = -1
    ix = 600
    plt.plot(delta/1e6, Mx[:,ix], label = 'Mx')
    plt.plot(delta/1e6, My[:,ix], label = 'My')
    plt.plot(delta/1e6, Mz[:,ix], label = 'Mz')
    plt.legend()
    plt.xlabel('Frequency (MHz)')

    plt.figure('Mx, My, Mz x spectrum')
    plt.plot(delta/1e6, Mx[:,ix]*spectrum.reshape(-1), label = 'Mx')
    plt.plot(delta/1e6, My[:,ix]*spectrum.reshape(-1), label = 'My')
    plt.plot(delta/1e6, Mz[:,ix]*spectrum.reshape(-1), label = 'Mz')
    plt.legend()
    plt.xlabel('Frequency (MHz)')

    plt.figure('Integrated frequency')
    plt.plot(times*1e6, np.sum(Mx*spectrum.reshape(-1,1), axis = 0), label = 'Mx')
    plt.plot(times*1e6, np.sum(My*spectrum.reshape(-1,1), axis = 0), label = 'My')
    plt.legend()
    plt.xlabel('Time (\u03bcs)')

    plt.figure('Integrated frequency Pure')
    plt.plot(times*1e6, np.sum(Mx,axis = 0), label = 'Mx')
    plt.plot(times*1e6, np.sum(My,axis = 0), label = 'My')
    plt.legend()
    plt.xlabel('Time (\u03bcs)')


    plt.figure('Spectrum')
    plt.plot(delta/1e6,spectrum)
    plt.xlabel('Frequency (MHz)')



    ## Plot results
    #plt.figure(figsize=(10, 6))
    #for i, d in enumerate(delta):
    #    plt.plot(times, results[i][2], label=f'? = {d:.2f}')
    #plt.xlabel('Time (s)')
    #plt.ylabel('Expectation value of Z')
    #plt.legend()
    #plt.title('Spin-1/2 Dynamics under Shaped Pulse with Relaxation')
    plt.show()
