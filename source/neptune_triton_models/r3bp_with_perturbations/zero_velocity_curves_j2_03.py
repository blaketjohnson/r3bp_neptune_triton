"""
zero_velocity_curves_j2_03.py

Plots zero-velocity curves (ZVCs) for the J2-perturbed CR3BP
model of the Neptune–Triton system.

Contours represent boundaries defined by the Jacobi constant,
marking energetically forbidden regions for a third body.

Author: Blake Johnson (based on Stuchi et al. 2008)

Usage:
    python zero_velocity_curves_j2_03.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# === Project Imports === #
from constants import J2_neptune, R_neptune_meters
from non_dimensionalizer import mu_nd, to_nondim_length
from lagrange_points_j2_02 import get_lagrange_points_and_constants

# === Parameters (nondimensional) === #
mu = mu_nd
R_nd = to_nondim_length(R_neptune_meters)
J2 = J2_neptune
J2_eff = J2 * R_nd**2
omega2 = 1 + 3 * J2_eff

# === Primaries === #
neptune = (-mu, 0)
triton = (1 - mu, 0)

# === Lagrange Points and Critical C Values === #
points = get_lagrange_points_and_constants()
L1 = points["L1"]
L2 = points["L2"]
L3 = points["L3"]
L4 = points["L4"]
L5 = points["L5"]
C_L1 = points["C_L1"]
C_L2 = points["C_L2"]
C_L3 = points["C_L3"]
C_L4 = points["C_L4"]
lagrange_points = {'L1': L1, 'L2': L2, 'L3': L3, 'L4': L4, 'L5': L5}

# === J2-Perturbed Effective Potential === #
def effective_potential(x, y, mu, J2_eff):
    r1 = np.sqrt((x + mu)**2 + y**2)
    r2 = np.sqrt((x - 1 + mu)**2 + y**2)
    U = (
        0.5 * omega2 * (x**2 + y**2)
        + (1 - mu) / r1
        + mu / r2
        + (1 - mu) * J2_eff / r1**3
    )
    return 2 * U  # Jacobi constant C = 2U

# === Mesh Grid === #
x = np.linspace(-1.5, 1.5, 800)
y = np.linspace(-1.2, 1.2, 600)
X, Y = np.meshgrid(x, y)
Z = effective_potential(X, Y, mu, J2_eff)

# === Plot Setup === #
combined_levels = sorted([C_L1, C_L2, C_L3, C_L4])
individual_Cs = [C_L1, C_L2, C_L3]
colors = ['red', 'green', 'blue', 'purple']
titles = [
    "Combined ZVCs",
    f"C = {C_L1:.6f} (L1 Open)",
    f"C = {C_L2:.6f} (L2 Open)",
    f"C = {C_L3:.6f} (L3 Open)"
]

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.flatten()

# === Top-Left Plot: Combined ZVCs === #
ax0 = axes[0]
contours = ax0.contour(X, Y, Z, levels=combined_levels, colors=colors, linewidths=1.5)
ax0.contourf(X, Y, Z, levels=[0, min(combined_levels)], colors=['lightgray'], alpha=0.5)
ax0.plot(*neptune, 'o', label='Neptune', color='blue')
ax0.plot(*triton, 'o', label='Triton', color='orange')
for label, (xpt, ypt) in lagrange_points.items():
    ax0.plot(xpt, ypt, 'k*', markersize=8)
    ax0.text(xpt + 0.03, ypt + 0.03, label, fontsize=8)
ax0.set_title(titles[0])
ax0.set_xlabel("x (nondimensional)")
ax0.set_ylabel("y (nondimensional)")
ax0.axis('equal')
ax0.grid(True)
legend_elements = [
    Line2D([0], [0], color='red', lw=2, label='C_L1'),
    Line2D([0], [0], color='green', lw=2, label='C_L2'),
    Line2D([0], [0], color='blue', lw=2, label='C_L3'),
    Line2D([0], [0], color='purple', lw=2, label='C_L4'),
    Line2D([0], [0], marker='o', color='w', label='Neptune', markerfacecolor='blue', markersize=8),
    Line2D([0], [0], marker='o', color='w', label='Triton', markerfacecolor='orange', markersize=8),
]
ax0.legend(handles=legend_elements, loc='lower left')

# === Remaining Subplots: Individual Levels === #
for i, (ax, C, title, color) in enumerate(zip(axes[1:], individual_Cs, titles[1:], colors)):
    cs = ax.contour(X, Y, Z, levels=[C], colors=color, linewidths=2.0)
    ax.contourf(X, Y, Z, levels=[0, C], colors=['lightgray'], alpha=0.5)
    ax.plot(*neptune, 'o', color='blue')
    ax.plot(*triton, 'o', color='orange')
    for label, (xpt, ypt) in lagrange_points.items():
        ax.plot(xpt, ypt, 'k*', markersize=8)
        ax.text(xpt + 0.03, ypt + 0.03, label, fontsize=8)
    ax.set_title(title)
    ax.set_xlabel("x (nondimensional)")
    ax.set_ylabel("y (nondimensional)")
    ax.axis('equal')
    ax.grid(True)

# === Final Layout === #
plt.tight_layout()
plt.savefig("zvc_2x2_combined_and_individual.png")
plt.show()









