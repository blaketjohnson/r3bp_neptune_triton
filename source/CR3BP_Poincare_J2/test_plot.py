import os
import numpy as np
import matplotlib.pyplot as plt

# === Parameters ===
CJ = 3.0148       # Jacobi constant
XI = -0.2         # x0 min
XF = 1.2          # x0 max
DX = 0.1000       # step between x0 values

# Folder with your generated .dat files
data_folder = "Poincare_data_CR3BP"

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

    # Match professor's filename format exactly
    fname = f"PY-C{CJ:.5f}Xi{round(x0, 5)}.dat"
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
plt.xlim(-1.1, 1.1)
plt.ylim(-1.2, 1.2)

# === Save plot ===
save_name = f"Poincare_C{CJ:.5f}_DX{DX:.4f}.png"
save_path = os.path.join(plot_folder, save_name)
plt.savefig(save_path, dpi=300)
print(f"✅ Plot saved to: {save_path}")

plt.show()

