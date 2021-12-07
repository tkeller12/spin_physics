import numpy as np
from matplotlib.pylab import *


I = 0.5 # spin quantum number
B0 = 20.
h_bar = 1.054571817e-34 #J / s
gamma = 42.577478518e6 # Hz / T
#T = np.r_[.001:300:100000j].reshape(-1, 1)
T = np.logspace(-3,2.5,100).reshape(-1, 1)
k_B = 1.38064852e-23

states = np.linspace(-1.*I, I, int(2.*I + 1)).reshape(1, -1)

energy = -1.*states*h_bar*gamma*B0


total_population = np.sum(np.exp(-energy / (k_B * T)), axis = 1).reshape(-1, 1)
#print(total_population)

populations = np.exp(-1.*energy / (k_B * T)) / total_population

polarization = 1 - populations

figure()
#plot(T, populations)
semilogx(T, populations)
#semilogx(T, polarization)
xlabel('Temperature (K)')
ylabel('P')
show()

