import numpy as np
from matplotlib.pylab import *

pts = 1024
omega = 100
tau = 0.2

t = np.r_[0:1:1j*pts]
t2 = np.r_[0:2:1j*2*pts]
t4 = np.r_[0:4:1j*4*pts]

y = np.exp(1j*omega*t) * np.exp(-t/tau)



y4 = np.hstack((y,y,y,y))
y2 = np.hstack((y,y))

ft = np.fft.fftshift(np.fft.fft(y))
ft2 = np.fft.fftshift(np.fft.fft(y2))
ft4 = np.fft.fftshift(np.fft.fft(y4))

dt = t[1] - t[0]
f = np.linspace(-0.5/dt,0.5/dt, num = pts, endpoint = False)
f2 = np.linspace(-0.5/dt,0.5/dt, num = 2*pts, endpoint = False)
f4 = np.linspace(-0.5/dt,0.5/dt, num = 4*pts, endpoint = False)

ft4 /= 4
ft2 /= 2

figure('time')
plot(t,y+2)

plot(t2,y2+1)
plot(t4,y4)

figure('ft')
plot(f,ft, '-')
plot(f2,ft2, '-')
plot(f4,ft4, '-')


show()
