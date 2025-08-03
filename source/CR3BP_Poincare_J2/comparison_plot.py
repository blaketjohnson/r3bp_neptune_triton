"""
plot_dual_professor_vs_yours.py

Panel 1: Professor style (low DX, discrete colors per x0 batch)
Panel 2: Your style (smooth gradient colormap)
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm

# === Parameters ===
CJ = 3.0148
XI = -0.2
XF = 1.2

DX_PROF = 0.1000  # from professor's plotting code
DX_YOUR = 0.005   # your normal resolution

# === Folders ===
prof_folder = "/Users/blakejohnson/Documents/r3bp_neptune_triton/CR3BP_Poincare/20231110P02"
your_folder = "Poincare_data"  # adjust if needed

# === Setup plotting ===
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['figure.autolayout'] = True

fig, axs = plt.subplots(1, 2, figsize=(14, 6), sharex=True, sharey=True)

# ------------------------------------------------
# PANEL 1 — Professor style with discrete colors
# ------------------------------------------------
prof_colors = plt.cm.tab20.colors  # 20 discrete colors
x0_values_prof = np.arange(XI, XF + DX_PROF / 2, DX_PROF)

for i, x0 in enumerate(x0_values_prof):
    # Match professor's exact filename pattern
    aux1 = str(round(CJ, 5))
    aux2 = str(round(x0, 5))
    fname = f"PY-C{aux1}Xi{aux2}.dat"
    path = os.path.join(prof_folder, fname)
    if not os.path.exists(path):
        continue

    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data[np.newaxis, :]

    # Cycle through discrete colors
    color = prof_colors[i % len(prof_colors)]
    axs[0].scatter(data[:, 0], data[:, 3], s=0.05, color=color)

axs[0].set_title("Professor Style", fontsize=16)
axs[0].set_xlabel(r"$x$ (ND)", fontsize=14)
axs[0].set_ylabel(r"$\dot{x}$ (ND)", fontsize=14)
axs[0].grid(True, alpha=0.3)

# ------------------------------------------------
# PANEL 2 — Your style (smooth gradient)
# ------------------------------------------------
cmap = cm.get_cmap("plasma")
norm = Normalize(vmin=XI, vmax=XF)

x0_values_your = np.arange(XI, XF + DX_YOUR / 2, DX_YOUR)

for x0 in x0_values_your:
    fname = f"PY-C{CJ:.5f}_Xi{x0:.5f}.dat"  # your naming pattern
    path = os.path.join(your_folder, fname)
    if not os.path.exists(path):
        continue

    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data[np.newaxis, :]

    axs[1].scatter(data[:, 0], data[:, 3],
                   s=0.2, alpha=0.5, color=cmap(norm(x0)))

axs[1].set_title("Your Style (J₂)", fontsize=16)
axs[1].set_xlabel(r"$x$ (ND)", fontsize=14)
axs[1].grid(True, alpha=0.3)

# Colorbar for your style
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=axs[1], label=r"$x_0$ (ND)")
cbar.ax.tick_params(labelsize=12)

# ------------------------------------------------
# Match axis limits
# ------------------------------------------------
for ax in axs:
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.2, 1.2)
    ax.ticklabel_format(style='plain', axis='x')
    ax.ticklabel_format(style='plain', axis='y')

# === Save & Show ===
os.makedirs("Plot_Results", exist_ok=True)
out_name = f"Comparison_C{CJ:.5f}.png"
save_path = os.path.join("Plot_Results", out_name)
plt.savefig(save_path, dpi=300)
print(f"✅ Dual plot saved to: {save_path}")

plt.show()
plt.close()