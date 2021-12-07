import numpy as np
from matplotlib.pylab import *

periods = 10
pts = 1024

#x = np.r_[0:periods*2*np.pi:1024j]
x = np.linspace(0, periods * 2 * np.pi, num = pts+1, endpoint = False)


omega = 5.0

y = np.cos(x*omega) - 1j* np.sin(x*omega)
y2 = np.hstack((y,y))

figure('stack')
plot(y2)
show()

ft = np.fft.fftshift(np.fft.fft(y))
ft2 = np.fft.fftshift(np.fft.fft(y2))

#data = np.zeros_like(ft)
#data[10] = 1
#data_ift = np.fft.ifft(data)

figure()
plot(x,y)

figure()
plot(np.real(ft))
plot(np.imag(ft))
show()

#figure()
#plot(x,np.real(data_ift), label = 'real')
#plot(x,np.imag(data_ift), label = 'imag')
#legend(loc = 'lower right')
#show()
