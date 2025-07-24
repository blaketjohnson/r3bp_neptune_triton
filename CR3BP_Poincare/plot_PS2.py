# Program to plot a Poincaré Section from file

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import cm
import numpy as np

# plt.rcParams["figure.figsize"] = [6.5, 6.5]
plt.rcParams["figure.autolayout"] = True
plt.rcParams['font.family'] = 'Times New Roman'
# plt.rcParams['font.sans-serif'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'

# plt.axvline(x = 4.45e-3, color = 'r')
# plt.axvline(x = (1 - 4.45e-3), color = 'b')

folder = "20231110P02"

CJ = 3.02

xi = 0.200
xf = 0.215
dx = 0.005
nx = round((xf - xi)/dx)

for ii in range(0, nx + 1):
    x0 = xi + float(ii)*dx

    aux1 = str(round(CJ, 5))
    aux2 = str(round(x0, 5))
    file = folder + '/' + 'PY-C' + aux1 + 'Xi' + aux2 + '.dat'

    ps = np.loadtxt(file, dtype=float)
    x = list(ps[:, 0])
    y = list(ps[:, 3])
    z = list(ps[:, 6])
    plt.scatter(x, y, s=0.05)

plt.xlabel(r'$x$ (nondimensional units)', fontsize=16)
plt.xticks(fontsize=14)
# plt.xlim(-10, 10)

plt.ylabel(r'$\dot{x}$ (nondimensional units)', fontsize=16)
plt.yticks(fontsize=14)

plt.show()