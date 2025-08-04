import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

# === Paths to the two files ===
prof_file = "/Users/blakejohnson/Documents/r3bp_neptune_triton/CR3BP_Poincare/20231110P02/PY-C3.0146Xi0.1.dat"

my_file = "/Users/blakejohnson/Documents/r3bp_neptune_triton/source/CR3BP_Poincare_J2/Poincare_data_CR3BP/PY-C3.01460Xi0.1.dat"


# === Load both ===
prof_data = np.loadtxt(prof_file)
my_data   = np.loadtxt(my_file)

if prof_data.ndim == 1:
    prof_data = prof_data[np.newaxis, :]
if my_data.ndim == 1:
    my_data = my_data[np.newaxis, :]

print(f"Professor crossings: {len(prof_data)}")
print(f"My code crossings:   {len(my_data)}")

# === Nearest-neighbor matching ===
tree = cKDTree(my_data[:, [0, 3]])  # x, xdot
dist, idx = tree.query(prof_data[:, [0, 3]], k=1)

matched_my = my_data[idx]

# === Differences ===
dx = prof_data[:, 0] - matched_my[:, 0]
dxdot = prof_data[:, 3] - matched_my[:, 3]

print(f"Mean Δx:     {np.mean(dx):.3e}")
print(f"Mean Δxdot:  {np.mean(dxdot):.3e}")
print(f"RMS Δx:      {np.sqrt(np.mean(dx**2)):.3e}")
print(f"RMS Δxdot:   {np.sqrt(np.mean(dxdot**2)):.3e}")
print(f"Max Δx:      {np.max(np.abs(dx)):.3e}")
print(f"Max Δxdot:   {np.max(np.abs(dxdot)):.3e}")

# === Overlay plot ===
plt.figure(figsize=(8,6))
plt.scatter(prof_data[:,0], prof_data[:,3], s=4, color='red', alpha=0.5, label="Professor")
plt.scatter(matched_my[:,0], matched_my[:,3], s=2, color='blue', alpha=0.5, label="Mine")
plt.xlabel(r"$x$ (ND)")
plt.ylabel(r"$\dot{x}$ (ND)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.title("Poincaré Crossing Comparison")
plt.show()

