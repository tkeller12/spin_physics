import numpy as np
from matplotlib.pylab import *

pts = 1024
omega = 100
tau = 0.1

t = np.r_[0:1:1j*pts]

y = np.exp(1j*omega*t) * np.exp(-t/tau)


#y = np.roll(y,100)



ft = np.fft.fftshift(np.fft.fft(y))

dt = t[1] - t[0]
f = np.linspace(-0.5/dt,0.5/dt, num = pts, endpoint = False)


figure('time')
plot(t,y)


figure('ft')
plot(f,ft, '-')


show()
