import numpy as np
from matplotlib.pylab import *

pts = 256
omega = 50
tau = 0.2

t = np.r_[0:1:1j*pts]
data = np.exp(1j*omega*t) * np.exp(-1*t/tau)

dt = t[1] - t[0]

f = np.linspace(-0.5/dt, 0.5/dt, pts, endpoint = False)

ft = np.fft.fftshift(np.fft.fft(data))


max_index = pts

fig, ax = plt.subplots(2,1)
#fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
init_index = 0
delta_index = 1
ax[0].plot(t, np.real(data))
ax[0].plot(t, np.imag(data))
l = ax[0].get_lines()
ax[1].plot(f, np.real(ft))
ax[1].plot(f, np.imag(ft))


l_ft = ax[1].get_lines()
ax[0].margins(x=0)

axcolor = "lightgoldenrodyellow"
axindex = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)

sindex = Slider(
    axindex,
    "index",
    -1 * max_index,
    max_index,
    valinit=init_index,
    valstep=delta_index,
)

def update(val):
    index = sindex.val
    updated_data = np.roll(data, index)
    updated_ft = np.fft.fftshift(np.fft.fft(updated_data))
    l[0].set_ydata(np.real(updated_data))
    l[1].set_ydata(np.imag(updated_data))
    l_ft[0].set_ydata(np.real(updated_ft))
    l_ft[1].set_ydata(np.imag(updated_ft))
    fig.canvas.draw_idle()

sindex.on_changed(update)

reset_ax = plt.axes([0.8, 0.025, 0.1, 0.04])
reset_button = Button(reset_ax, "Reset", color=axcolor, hovercolor="0.975")

inc_ax = plt.axes([0.6, 0.025, 0.1, 0.04])
inc_button = Button(inc_ax, "+", color=axcolor, hovercolor="0.975")

dec_ax = plt.axes([0.4, 0.025, 0.1, 0.04])
dec_button = Button(dec_ax, "-", color=axcolor, hovercolor="0.975")

def reset(event):
    sindex.reset()

def inc(event):
    sindex.set_val(sindex.val + 1)

def dec(event):
    sindex.set_val(sindex.val - 1)

reset_button.on_clicked(reset)
inc_button.on_clicked(inc)
dec_button.on_clicked(dec)

plt.show()
index = sindex.val

