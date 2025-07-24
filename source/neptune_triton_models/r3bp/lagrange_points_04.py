"""
lagrange_points_04.py

This script computes the five Lagrange points (L1 through L5) in the Circular Restricted Three-Body Problem (CR3BP)
for the Neptune–Triton system. These are positions in a rotating frame where gravitational and centrifugal forces
balance.

References:
- Murray & Dermott, Solar System Dynamics, 1999, Ch. 3
- Gradient of the effective potential: ∇U = 0 (Murray p. 67, Eq. 3.23)
- Lagrange point locations: solutions to ∇U = 0

Usage:
    python lagrange_points_04.py
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root
import sys

# Import necessary parameters and functions
sys.path.append('/Users/blakejohnson/Documents/rtbp_neptune_triton/neptune_triton_models/r3bp')
from r3bp_calculations_02 import mu_r3bp, gradient_U  # updated for modular naming

def find_Lagrange_points(mu):
    """
    Computes all five Lagrange points (L1–L5) for a given mass ratio mu.

    Parameters:
        mu : float
            Mass parameter, defined as mu = m2 / (m1 + m2)

    Returns:
        tuple of np.array
            The coordinates of L1, L2, L3, L4, and L5 in rotating frame (x, y, z)
    """

    def find_L1(mu):
        guess = [0.5 - mu, 0.0, 0.0]
        sol = root(lambda coords: gradient_U(*coords, mu), guess)
        if sol.success:
            return sol.x
        raise RuntimeError(f"Failed to find L1: {sol.message}")

    def find_L2(mu):
        guess = [1.5 - mu, 0.0, 0.0]
        sol = root(lambda coords: gradient_U(*coords, mu), guess)
        if sol.success:
            return sol.x
        raise RuntimeError(f"Failed to find L2: {sol.message}")

    def find_L3(mu):
        guess = [-1.0 - mu, 0.0, 0.0]
        sol = root(lambda coords: gradient_U(*coords, mu), guess)
        if sol.success:
            return sol.x
        raise RuntimeError(f"Failed to find L3: {sol.message}")

    def find_L4(mu):
        guess = [0.5 - mu, np.sqrt(3)/2, 0.0]
        sol = root(lambda coords: gradient_U(*coords, mu), guess)
        if sol.success:
            return sol.x
        raise RuntimeError(f"Failed to find L4: {sol.message}")

    def find_L5(mu):
        guess = [0.5 - mu, -np.sqrt(3)/2, 0.0]
        sol = root(lambda coords: gradient_U(*coords, mu), guess)
        if sol.success:
            return sol.x
        raise RuntimeError(f"Failed to find L5: {sol.message}")

    L1 = find_L1(mu)
    L2 = find_L2(mu)
    L3 = find_L3(mu)
    L4 = find_L4(mu)
    L5 = find_L5(mu)
    print(f"Lagrange Points (CR3BP):\nL1: {L1}\nL2: {L2}\nL3: {L3}\nL4: {L4}\nL5: {L5}")
    return L1, L2, L3, L4, L5

def plot_lagrange_points(L_points):
    """
    Plot Neptune, Triton, and Lagrange points in the rotating frame.
    :param L_points: tuple of np.array, the Lagrange points (L1, L2, L3, L4, L5)
    """
    plt.figure(figsize=(10, 8))

    # Primary and secondary bodies in rotating frame
    neptune = np.array([0.0, 0.0])
    triton = np.array([1.0, 0.0])

    plt.scatter(*neptune, color='blue', label='Neptune')
    plt.scatter(*triton, color='saddlebrown', label='Triton')

    for i, point in enumerate(L_points, start=1):
        color = 'red' if i == 4 else 'black'
        plt.scatter(point[0], point[1], color=color, label=f'L{i}')
        plt.text(point[0], point[1], f'L{i}', fontsize=12, ha='center', va='bottom')

    plt.text(*neptune, ' Neptune', fontsize=12, ha='center', va='top')
    plt.text(*triton, ' Triton', fontsize=12, ha='center', va='top')

    plt.xlabel('X [unitless]')
    plt.ylabel('Y [unitless]')
    plt.title('Lagrange Points in the Neptune–Triton System (CR3BP)')
    plt.grid(True)
    plt.legend()
    plt.axis('equal')
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    # Compute Lagrange points
    L1, L2, L3, L4, L5 = find_Lagrange_points(mu_r3bp)
    plot_lagrange_points((L1, L2, L3, L4, L5))





