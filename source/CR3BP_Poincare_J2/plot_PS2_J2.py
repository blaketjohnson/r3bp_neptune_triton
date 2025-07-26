"""
plot_single_Poincare_J2.py

Plots a single Poincaré section at a chosen Jacobi constant (CJ)
and across a small set of x0 values. Designed for inspecting local structure
without rerunning long simulations. Results are saved to Plot_Results.

Author: Blake T. Johnson
Summer Thesis Project, Phase 1
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm

# === Plot Target Parameters ===
CJ = 2.990                  # Jacobi constant to plot
XI = 0.200                  # Starting x0 (ND)
XF = 0.215                  # Ending x0 (ND)
DX = 0.001                  # Step size for x0 (ND)

# === Paths and Style ===
data_folder = "Poincare_data"
plot_folder = "Plot_Results"
os.makedirs(plot_folder, exist_ok=True)

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['figure.autolayout'] = True

# === Derived Range ===
x0_values = np.arange(XI, XF + DX/2, DX)
fig, ax = plt.subplots(figsize=(8, 6))
cmap = cm.get_cmap("plasma")
norm = Normalize(vmin=XI, vmax=XF)

# === Plot Loop ===
for x0 in x0_values:
    fname = f"PY-C{CJ:.5f}_Xi{x0:.5f}.dat"
    path = os.path.join(data_folder, fname)
    if not os.path.exists(path):
        continue
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    ax.scatter(data[:, 0], data[:, 3], s=0.5, color=cmap(norm(x0)))

# === Labels and Save ===
ax.set_xlabel(r'$x$ (ND)', fontsize=14)
ax.set_ylabel(r'$\dot{x}$ (ND)', fontsize=14)
ax.set_title(f"Poincaré Section at C = {CJ:.5f}", fontsize=16)
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
plt.colorbar(sm, ax=ax, label=r'$x_0$ (ND)')
ax.grid(True)

out_name = f"Poincare_C{CJ:.5f}.png"
plt.savefig(os.path.join(plot_folder, out_name), dpi=300)
plt.close()



