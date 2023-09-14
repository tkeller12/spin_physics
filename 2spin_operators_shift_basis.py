import numpy as np
from matplotlib.pylab import *
from matplotlib import cm

sigma_x = 0.5*np.r_[[[0, 1],[1, 0]]]
sigma_y = 0.5*np.r_[[[0,-1j],[1j, 0]]]
sigma_z = 0.5*np.r_[[[1, 0],[0, -1]]]
Identity = np.eye(2)


Sx = np.kron(sigma_x, Identity)
Sy = np.kron(sigma_y, Identity)
Sz = np.kron(sigma_z, Identity)

Ix = np.kron(Identity, sigma_x)
Iy = np.kron(Identity, sigma_y)
Iz = np.kron(Identity, sigma_z)

SxIz = np.dot(Sx, Iz)
SyIz = np.dot(Sy, Iz)
SzIz = np.dot(Sz, Iz)

SzIx = np.dot(Sz, Ix)
SzIy = np.dot(Sz, Iy)

SxIx = np.dot(Sx, Ix)
SyIy = np.dot(Sx, Ix)

SxIy = np.dot(Sx, Iy)
SyIx = np.dot(Sy, Ix)

operators = {
        'E': np.eye(4),
        '$S_0$': np.sqrt(2)*Sz,
        '$I_0$': np.sqrt(2)*Iz,
        '$S_0I_0$' : SzIz,

        '$S^+$': Sx + 1j*Sy,
        '$S^-$': Sx - 1j*Sy,
        '$S^+I_0$' : np.dot(Sx + 1j*Sy,Iz),
        '$S^-I_0$' : np.dot(Sx - 1j*Sy,Iz),

        '$I^+$': Ix + 1j*Iy,
        '$I^-$': Ix - 1j*Iy,
        '$S_0I^+$' : np.dot(Sz,Ix + 1j*Iy),
        '$S_0I^-$' : np.dot(Sz,Ix - 1j*Iy),

        '$S^+I^+$' : np.dot(Ix + 1j*Iy,Sx + 1j*Sy),
        '$S^-I^+$' : np.dot(Ix - 1j*Iy,Sx + 1j*Sy),
        '$S^+I^-$' : np.dot(Ix + 1j*Iy,Sx - 1j*Sy),
        '$S^-I^-$' : np.dot(Ix - 1j*Iy,Sx - 1j*Sy),
        }

fig, axs = plt.subplots(4, 4)

ix = 0

ops = list(operators.keys())
for row in range(4):
    for column in range(4):
        O = ops[ix]
        axs[row, column].matshow(np.abs(operators[O]), cmap = cm.YlGn)
#        axs[row, column].matshow(np.imag(operators[O]), cmap = cm.YlGn)
        axs[row, column].set_title(O)
        axs[row, column].get_xaxis().set_ticks([])
        axs[row, column].get_yaxis().set_ticks([])
        ix += 1

tight_layout()
show()
