# Program to plot a Poincaré Section from file

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import cm
import numpy as np
from CR3BP_Poincare_j2 import initial_parameters
import os
from numba import njit

# plt.rcParams["figure.figsize"] = [6.5, 6.5]
plt.rcParams["figure.autolayout"] = True
plt.rcParams['font.family'] = 'Times New Roman'
# plt.rcParams['font.sans-serif'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'

# plt.axvline(x = 4.45e-3, color = 'r')
# plt.axvline(x = (1 - 4.45e-3), color = 'b')

#### Reading the input parameters from the class ####
init_par = initial_parameters()

folder = init_par.folder

C_range = init_par.C_range
C0 = float(C_range[0])
CF = float(C_range[1])
dC = float(C_range[2])

CJ = C0

x_range = init_par.x_range
xi = x_range[0]
xf = x_range[1]
dx = x_range[2]

nx = round((xf - xi)/dx)

for ii in range(0, nx + 1):
    x0 = xi + float(ii)*dx

    aux1 = "%05.3f" % CJ
    aux2 = "%05.3f" % x0

    # Test if the file exists
    file = folder + '/' + 'PY-C' + aux1 + 'Xi' + aux2 + '.dat'
    if not os.path.isfile(file):
        print(f"File {file} does not exist. Skipping...")
        continue

    ps = np.loadtxt(file, dtype=float)
    x = list(ps[:, 0])
    y = list(ps[:, 3])
    #z = list(ps[:, 6])
    plt.scatter(x, y, s=0.05) #, c='blue'
    # plt.scatter(x, y, s=0.05, alpha=0.7, c=z, cmap='viridis', norm=mcolors.Normalize(vmin=-1.0, vmax=1.0))

plt.xlabel(r'$x$ (nondimensional units)', fontsize=16)
plt.xticks(fontsize=14)

plt.title(r'Poincaré Section for $C = $' + "%05.3f" % CJ, fontsize=16)
# plt.xlim(-1.0, 1.0)
# plt.ylim(-1.0, 1.0)

# Set plot size and aspect ratio
# plt.gcf().set_size_inches(6.5, 6.5)
# plt.gca().set_aspect('equal', adjustable='box')
# plt.gca().set_aspect('auto', adjustable='box')
# plt.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
plt.gca().set_aspect('auto', adjustable='box')

# plt.xlim(-10, 10)

plt.ylabel(r'$\dot{x}$ (nondimensional units)', fontsize=16)
plt.yticks(fontsize=14)

plt.show()