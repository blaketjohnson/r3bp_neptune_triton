"""
plot_single_Poincare_J2_matched.py

Plots a single Poincaré section for a chosen Jacobi constant (CJ)
and x0 sweep range, with style/formatting matched closely to professor's plots.

Author: Blake T. Johnson
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm

# === Plot Target Parameters ===
CJ = 2.9000       # Jacobi constant to plot
XI = -0.2         # Start x0 (ND)
XF = 1.2          # End x0 (ND)
DX = 0.005        # Step size for x0 (ND)

# === Paths ===
data_folder = "Poincare_data"
plot_folder = "Plot_Results"
os.makedirs(plot_folder, exist_ok=True)

# === Style ===
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['figure.autolayout'] = True

# === Prepare figure ===
fig, ax = plt.subplots(figsize=(8, 6))
cmap = cm.get_cmap("plasma")  # or "viridis" if you want cooler tones
norm = Normalize(vmin=XI, vmax=XF)

# === Sweep range ===
x0_values = np.arange(XI, XF + DX/2, DX)
colors_used = []

# === Plot loop ===
for x0 in x0_values:
    fname = f"PY-C{CJ:.5f}_Xi{x0:.5f}.dat"
    path = os.path.join(data_folder, fname)
    if not os.path.exists(path):
        continue

    data = np.loadtxt(path)
    if data.ndim == 1:  # Ensure 2D array
        data = data[np.newaxis, :]

    # Filter absurd velocities (ND sanity check)
    mask = np.abs(data[:, 3]) < 10.0
    data = data[mask]
    if data.size == 0:
        continue

    colors_used.append(x0)
    ax.scatter(data[:, 0], data[:, 3],
           s=0.2, alpha=0.5, color=cmap(norm(x0)))  # visibility tuned for dense regions

# === Colorbar ===
if colors_used:
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, label=r'$x_0$ (ND)')
    cbar.ax.tick_params(labelsize=12)

# === Labels & Title ===
ax.set_xlabel(r'$x$ (ND)', fontsize=14)
ax.set_ylabel(r'$\dot{x}$ (ND)', fontsize=14)
ax.set_title(f"Poincaré Section at C = {CJ:.5f}", fontsize=16)
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)


# Match professor’s clean axis look
ax.ticklabel_format(style='plain', axis='x')
ax.ticklabel_format(style='plain', axis='y')
ax.grid(True, alpha=0.3)

# === Save & Show ===
out_name = f"Poincare_C{CJ:.5f}.png"
save_path = os.path.join(plot_folder, out_name)
plt.savefig(save_path, dpi=300)
print(f"✅ Plot saved to: {save_path}")

plt.show()
plt.close()






