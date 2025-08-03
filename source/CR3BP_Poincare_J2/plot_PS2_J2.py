"""
plot_single_Poincare_J2_filtered.py

Same as your Poincaré plotting script but skips points with absurd ND velocities
to prevent them from stretching the y-axis.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm

# === Plot Target Parameters ===
CJ = 3.0148
XI = -0.2
XF = 1.2
DX = 0.1000

# === Paths and Style ===
data_folder = "Poincare_data"
plot_folder = "Plot_Results"
os.makedirs(plot_folder, exist_ok=True)

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['figure.autolayout'] = True

# === Prepare figure and colormap ===
x0_values = np.arange(XI, XF + DX/2, DX)
fig, ax = plt.subplots(figsize=(8, 6))
cmap = cm.get_cmap("viridis")

colors_used = []

# === Plot Loop with filtering ===
for x0 in x0_values:
    fname = f"PY-C{CJ:.5f}_Xi{x0:.5f}.dat"
    path = os.path.join(data_folder, fname)
    if not os.path.exists(path):
        continue

    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data[np.newaxis, :]

    # Debug print of ND velocity range
    print(f"x0 = {x0:.3f} -> vx min/max = {data[:,3].min():.4f} / {data[:,3].max():.4f}")

    # === Velocity filter: remove bad points ===
    mask = np.abs(data[:, 3]) < 10.0  # ND velocities should be small (~<2)
    data = data[mask]
    if data.size == 0:
        continue

    colors_used.append(x0)
    ax.scatter(data[:, 0], data[:, 3],
               s=3, alpha=0.8,
               color=cmap((x0 - XI) / (XF - XI)))

# === Colorbar ===
if colors_used:
    norm = Normalize(vmin=min(colors_used), vmax=max(colors_used))
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, label=r'$x_0$ (ND)')

# === Labels and formatting ===
ax.set_xlabel(r'$x$ (ND)', fontsize=14)
ax.set_ylabel(r'$\dot{x}$ (ND)', fontsize=14)
ax.set_title(f"Poincaré Section at C = {CJ:.5f}", fontsize=16)
ax.grid(True)
ax.ticklabel_format(style='plain', axis='x')
ax.ticklabel_format(style='plain', axis='y')

# === Save and show ===
out_name = f"Poincare_C{CJ:.5f}_filtered.png"
save_path = os.path.join(plot_folder, out_name)
plt.savefig(save_path, dpi=300)
print(f"\n✅ Poincaré section saved to: {save_path}")

plt.show()
plt.close()





