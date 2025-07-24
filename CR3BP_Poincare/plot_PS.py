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

# cm = plt.cm.get_cmap('viridis')

ps = np.loadtxt("20231110P02/PY-C3.02Xi0.2.dat", dtype=float)
x = list(ps[:, 0])
y = list(ps[:, 3])
z = list(ps[:, 6])

sc = plt.scatter(x, y, c=z, vmin=0.8, vmax=1.2, s=0.2, cmap=cm)
plt.colorbar(sc)

plt.xlabel(r'$x$ (nondimensional units)', fontsize=16)
plt.xticks(fontsize=14)
# plt.xlim(-10, 10)

plt.ylabel(r'$\dot{x}$ (nondimensional units)', fontsize=16)
plt.yticks(fontsize=14)


plt.show()