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
import glob
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from parameters import *  # Pulls mapping_mode, J2_enabled, DX, XI, XF, etc.

# --- Determine CJ from parameters ---
if C0 == CF:
    CJ = C0
else:
    CJ = C0  # Default to first value in sweep

filter_by_x_range = True  # Only include files within XI → XF

# --- Timestamp for saved plot ---
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

# --- Base project path ---
project_root = "/Users/blakejohnson/Documents/r3bp_neptune_triton/source/CR3BP_Poincare_J2"

# --- Build folder paths ---
mode_folder = "highres" if mapping_mode.lower() == "highres" else "global"
perturb_folder = "perturbed" if J2_enabled else "non_perturbed"

data_folder = os.path.join(project_root, mode_folder, perturb_folder, "data")
results_folder = os.path.join(project_root, mode_folder, perturb_folder, "results")

os.makedirs(results_folder, exist_ok=True)

print(f"Looking for data in: {data_folder}")
print(f"Plots will be saved in: {results_folder}")

# --- Search for matching .dat files ---
pattern = f"CJ{CJ:.5f}_Xi*_DX{DX:.4f}_*_{mapping_mode.lower()}_*.dat"
files = sorted(glob.glob(os.path.join(data_folder, pattern)))

if not files:
    print("No matching .dat files found.")
    exit()

# --- Optional: Filter by x-range ---
if filter_by_x_range:
    filtered_files = []
    for f in files:
        try:
            # Extract Xi value from filename (format: ..._Xi<value>_...)
            xi_str = [part for part in os.path.basename(f).split("_") if part.startswith("Xi")][0]
            xi_val = float(xi_str[2:])
            if XI <= xi_val <= XF:
                filtered_files.append(f)
        except (IndexError, ValueError):
            continue
    files = filtered_files

if not files:
    print("No files match the specified x-range.")
    exit()

# --- Matplotlib style ---
plt.rcParams["figure.autolayout"] = True
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "cm"

# --- Prepare arrays for all points ---
all_x = []
all_xdot = []
all_xi_vals = []

for path in files:
    data = np.loadtxt(path, dtype=float)
    if data.size == 0:
        continue
    if data.ndim == 1:
        data = data[np.newaxis, :]  # ensure 2D

    # Extract initial x0 from filename for coloring
    try:
        xi_str = [part for part in os.path.basename(path).split("_") if part.startswith("Xi")][0]
        xi_val = float(xi_str[2:])
    except (IndexError, ValueError):
        xi_val = np.nan

    # Append all points with their corresponding xi_val
    all_x.extend(data[:, 0])
    all_xdot.extend(data[:, 3])
    all_xi_vals.extend([xi_val] * len(data))

# Convert to numpy arrays
all_x = np.array(all_x)
all_xdot = np.array(all_xdot)
all_xi_vals = np.array(all_xi_vals)

# --- Single scatter plot for consistent color mapping ---
sc = plt.scatter(all_x, all_xdot, s=0.05, c=all_xi_vals,
                 cmap="plasma", vmin=XI, vmax=XF)

# --- Labels, title, and colorbar ---
title_str = (
    f"Poincaré Map — CJ={CJ:.5f}, DX={DX:.4f}, "
    f"x_0 Range: {XI:.3f} → {XF:.3f}\n"
    f"Mode={mapping_mode.capitalize()}, J_2={'On' if J2_enabled else 'Off'}"
)
plt.title(title_str, fontsize=14)

cbar = plt.colorbar(sc)
cbar.set_label(r"Initial $x_0$ (nondimensional)", fontsize=12)

plt.xlabel(r"$x$ (nondimensional units)", fontsize=16)
plt.ylabel(r"$\dot{x}$ (nondimensional units)", fontsize=16)
plt.grid(True, alpha=0.3)

# --- Adjust margins ---
x_margin = 0.05 * (np.max(all_x) - np.min(all_x))
y_margin = 0.05 * (np.max(all_xdot) - np.min(all_xdot))
plt.xlim(np.min(all_x) - x_margin, np.max(all_x) + x_margin)
plt.ylim(np.min(all_xdot) - y_margin, np.max(all_xdot) + y_margin)

# --- Save plot ---
# === Perturbation label for filename ===
perturb_tag = "perturbed" if J2_enabled else "non_perturbed"

# === New filename format ===
save_name = f"C{CJ:.5f}_x0_{XI:.4f}_xf_{XF:.4f}_{perturb_tag}.png"
save_path = os.path.join(results_folder, save_name)

plt.savefig(save_path, dpi=300)

print(f"Plot saved to: {save_path}")
print(f"Total files plotted: {len(files)}")
plt.show()







