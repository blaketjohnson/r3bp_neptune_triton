"""
zero_velocity_curves_05.py

This script computes and visualizes Zero-Velocity Curves (ZVCs) in the Circular 
Restricted Three-Body Problem (CR3BP) for the Neptune–Triton system. These curves 
represent boundaries in the rotating frame where a test particle's velocity would be zero 
for a given Jacobi constant, thus defining regions forbidden to motion.

References:
- Murray & Dermott, Solar System Dynamics, Cambridge University Press, 1999
- Jacobi Integral: (Murray p. 68, Eq. 3.23)
- Zero-Velocity Curves: (Murray p. 102, Sec. 3.10)

Usage:
    python zero_velocity_curves_05.py

Functionality:
- Computes Jacobi constants for the collinear Lagrange points (L1, L2, L3)
- Evaluates and plots zero-velocity curves (ZVCs) in the x-y plane for each C value
- Marks Neptune, Triton, and all five Lagrange points on the plot

Assumptions:
- Motion is restricted to the plane (z = 0)
- Mass ratio (mu) is based on the Neptune–Triton system
- Unitless distance and time normalized by Neptune–Triton parameters
"""


import numpy as np
from scipy.optimize import root
import matplotlib.pyplot as plt
from lagrange_points_04 import find_Lagrange_points
import sys

# Import necessary parameters and functions
sys.path.append('/Users/blakejohnson/Documents/rtbp_neptune_triton/neptune_triton_models/r3bp')
from r3bp_calculations_02 import mu_r3bp, gradient_U  # updated for modular naming


def U(x, y, z, mu):
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x + mu - 1)**2 + y**2 + z**2)
    return 0.5 * (x**2 + y**2 + z**2) + (1 - mu) / r1 + mu / r2

def gradient_U(x, y, z, mu):
    r1 = (x + mu)**2 + y**2 + z**2
    r2 = (x + mu - 1)**2 + y**2 + z**2
    
    Ux = x - (1 - mu) * (x + mu) / r1**1.5 - mu * (x + mu - 1) / r2**1.5
    Uy = y - (1 - mu) * y / r1**1.5 - mu * y / r2**1.5
    Uz = -(1 - mu) * z / r1**1.5 - mu * z / r2**1.5
    
    return np.array([Ux, Uy, Uz])

def jacobi_constant(x, y, z, mu):
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x + mu - 1)**2 + y**2 + z**2)
    return x**2 + y**2 + 2 * (1 - mu) / r1 + 2 * mu / r2

def zero_velocity_curve(x, y, C, mu):
    z = 0  # Zero-velocity curves in the x-y plane
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x + mu - 1)**2 + y**2 + z**2)
    return x**2 + y**2 + 2 * (1 - mu) / r1 + 2 * mu / r2 - C

# Calculate Lagrange points
L1, L2, L3, L4, L5 = find_Lagrange_points(mu_r3bp)

# Calculate Jacobi constants for L1, L2, and L3
C_L1 = jacobi_constant(L1[0], L1[1], L1[2], mu_r3bp)
C_L2 = jacobi_constant(L2[0], L2[1], L2[2], mu_r3bp)
C_L3 = jacobi_constant(L3[0], L3[1], L3[2], mu_r3bp)

# Grid for plotting zero-velocity curves
x_vals = np.linspace(-2, 2, 400)
y_vals = np.linspace(-2, 2, 400)
X, Y = np.meshgrid(x_vals, y_vals)

Z_L1 = zero_velocity_curve(X, Y, C_L1, mu_r3bp)
Z_L2 = zero_velocity_curve(X, Y, C_L2, mu_r3bp)
Z_L3 = zero_velocity_curve(X, Y, C_L3, mu_r3bp)

# Plot Neptune, Triton, Lagrange points, and zero-velocity curves in 2D
neptune_position = np.array([0.0, 0.0])
triton_position = np.array([1.0, 0.0])

if __name__ == '__main__':
    plt.figure(figsize=(10, 8))
    plt.scatter(neptune_position[0], neptune_position[1], color='b', label='Neptune')
    plt.scatter(triton_position[0], triton_position[1], color='g', label='Triton')
    plt.scatter(L1[0], L1[1], color='r', label='L1')
    plt.scatter(L2[0], L2[1], color='orange', label='L2')
    plt.scatter(L3[0], L3[1], color='purple', label='L3')
    plt.scatter(L4[0], L4[1], color='cyan', label='L4')
    plt.scatter(L5[0], L5[1], color='magenta', label='L5')

    plt.text(neptune_position[0], neptune_position[1], ' Neptune', fontsize=12, ha='center', va='bottom')
    plt.text(triton_position[0], triton_position[1], ' Triton', fontsize=12, ha='center', va='bottom')
    plt.text(L1[0], L1[1], ' L1', fontsize=12, ha='center', va='bottom')
    plt.text(L2[0], L2[1], ' L2', fontsize=12, ha='center', va='bottom')
    plt.text(L3[0], L3[1], ' L3', fontsize=12, ha='center', va='bottom')
    plt.text(L4[0], L4[1], ' L4', fontsize=12, ha='center', va='bottom')
    plt.text(L5[0], L5[1], ' L5', fontsize=12, ha='center', va='bottom')
    plt.text(0, 0, 'O', fontsize=12, ha='right', va='top')  # Origin (center of mass)

    plt.contour(X, Y, Z_L1, levels=[0], colors='r', linestyles='--', linewidths=0.5)
    plt.contour(X, Y, Z_L2, levels=[0], colors='orange', linestyles='--', linewidths=0.5)
    plt.contour(X, Y, Z_L3, levels=[0], colors='purple', linestyles='--', linewidths=0.5)

    plt.xlabel('Distance in unitless dimensions')
    plt.ylabel('Distance in unitless dimensions')
    plt.title('Positions of Lagrange Points and Zero-Velocity Curves with Neptune and Triton')
    plt.legend()
    plt.grid(True)

    # Adjust the limits to zoom out
    plt.xlim(-1.1, 1.1)
    plt.ylim(-1.1, 1.1)

    plt.show()

