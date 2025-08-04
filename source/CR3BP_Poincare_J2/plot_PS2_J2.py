"""
plot_CR3BP_mapper.py
---------------------------------------
Plots Poincaré surface-of-section maps from pre-generated CR3BP
data files (.dat) produced by CR3BP_Poincare_J2.py.

Features:
    • User-selectable mapping mode ("global" or "highres") at runtime
    • Reads data from mode-specific folders:
        Poincare_data_global/  or  Poincare_data_highres/
    • Matches solver’s .dat naming format:
        PY-C<Jacobi>Xi<Initial_x0>_<mode>.dat
    • Saves plot images with mode tag in filename, e.g.:
        Poincare_C3.00000_DX0.1000_global.png
    • Preserves professional plot styling (Times New Roman, LaTeX math formatting)

Usage:
    1. Set the plotting parameters in this script:
        - CJ : Jacobi constant to plot
        - XI, XF : x₀ sweep range
        - DX : step size
        - mapping_mode : "global" or "highres"
    2. Ensure .dat files exist for the chosen parameters/mode.
    3. Run:
        $ python plot_CR3BP_mapper.py
    4. The script will:
        - Read matching .dat files from the mode’s folder
        - Plot the Poincaré section
        - Save the figure to:
            Plot_Results/Poincare_C<Jacobi>_DX<step>_<mode>.png

Author:
    Blake T. Johnson
    Thesis Project
    (c) 2025
"""


import os
import numpy as np
import matplotlib.pyplot as plt

# === User selects which data to plot ===
CJ = 2.9000        # Jacobi constant to plot
XI = -0.2          # x0 min
XF = 1.2           # x0 max
DX = 0.1000        # step between x0 values
# === Mapping Mode ===
# "highres" = event detection, fine DX, exact crossings
# "global"  = fixed step, midpoint interpolation, professor-style
mapping_mode = "global"  

# Folder with your generated .dat files
data_folder = f"Poincare_data_{mapping_mode.lower()}"


# Folder to save plot
plot_folder = "Plot_Results"
os.makedirs(plot_folder, exist_ok=True)

# === Matplotlib style ===
plt.rcParams["figure.autolayout"] = True
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "cm"

# === Compute number of steps ===
nx = round((XF - XI) / DX)

# === Plot loop ===
for ii in range(0, nx + 1):
    x0 = XI + float(ii) * DX

    fname = f"PY-C{CJ:.5f}Xi{round(x0, 5)}_{mapping_mode.lower()}.dat"
    path = os.path.join(data_folder, fname)

    if not os.path.isfile(path):
        continue

    data = np.loadtxt(path, dtype=float)

    # Skip empty or bad files
    if data.size == 0:
        continue
    if data.ndim == 1:
        data = data[np.newaxis, :]

    # Plot x vs xdot (cols 0 and 3)
    plt.scatter(data[:, 0], data[:, 3], s=0.05)  # default discrete color cycle

# === Labels ===
plt.xlabel(r"$x$ (nondimensional units)", fontsize=16)
plt.ylabel(r"$\dot{x}$ (nondimensional units)", fontsize=16)
plt.grid(True, alpha=0.3)

# === Match professor's zoom ===
plt.xlim(-1.25, 1.25)
plt.ylim(-1.25, 1.25)

# === Save plot ===
save_name = f"Poincare_C{CJ:.5f}_DX{DX:.4f}_{mapping_mode.lower()}.png"
save_path = os.path.join(plot_folder, save_name)
plt.savefig(save_path, dpi=300)
print(f"Plot saved to: {save_path}")

plt.show()






