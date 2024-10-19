import numpy as np
from scipy.linalg import expm
from matplotlib.pylab import *
from pyDEER import awg


resolution = 1e-9

tp = 500e-9
t,pulse = awg.plane_wave(tp,50e6, resolution = resolution)
#pulse += 0.2 + 0.5j

t, sinc_pulse = awg.sinc(tp, 5, resolution = resolution)

pulse *= sinc_pulse


wave = np.hstack((np.zeros(1000),pulse,np.zeros(1000)))

pts = len(wave)

ft_wave = np.fft.fftshift(np.fft.fft(wave))

#filter_bandwidth = [3e6, 5e6, 10e6, 20e6, 50e6]
filter_bandwidth = [3e6, 10e6, 50e6]


t = np.linspace(0,pts*resolution,num = pts)

f = np.linspace(-0.5/resolution, 0.5/resolution, num = pts)

figure()
for ix,bw in enumerate(filter_bandwidth):
    filter_points = np.abs(f) < (bw/2)
    temp_ft_wave = ft_wave

    temp_ft_wave[filter_points] = 0


    filtered_wave = np.fft.ifft(np.fft.fftshift(temp_ft_wave))

    title('Simulated Pulse Shape')
    plot(t/1e-9, np.real(filtered_wave),label = '%i MHz BW'%(bw/1e6))
    legend()
    xlabel('Time (ns)')
    ylabel('Signal (a.u.)')
    xlim(np.min(t/1e-9),np.max(t/1e-9))
#    plot(t, np.imag(filtered_wave))
show()

#figure()
#plot(t, wave)
#
#figure()
#plot(f, ft_wave)

