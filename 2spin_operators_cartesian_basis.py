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
        '$S_z$': Sz,
        '$I_z$': Iz,
        '$S_zI_z$' : SzIz,

        '$S_x$': Sx,
        '$S_y$': Sy,
        '$S_xI_z$' : SxIz,
        '$S_yI_z$' : SyIz,

        '$S_zI_x$' : SzIx,
        '$S_zI_y$' : SzIy,
        '$I_x$': Ix,
        '$I_y$': Iy,


        '$S_xI_x$' : SxIx,
        '$S_yI_y$' : SyIy,
        '$S_xI_y$' : SxIy,
        '$S_yI_x$' : SyIx,
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
