"""
stability_lagrange_j2_11.py

This script analyzes the stability of the Lagrange Points (L1 through L5) in the 
J2-perturbed Circular Restricted Three-Body Problem (CR3BP) for the Neptune–Triton system.

Stability is assessed using eigenvalue analysis of the linearized equations near each point.

Usage:
    python stability_lagrange_j2_05.py
"""

import numpy as np
import sys

from scipy.linalg import eigvals
from constants import J2_neptune, R_neptune_meters, a_triton_meters
from lagrange_points_j2_02 import get_lagrange_points_and_constants

sys.path.append('/Users/blakejohnson/Documents/r3bp_neptune_triton/source/neptune_triton_models/r3bp')
from r3bp_calculations_02 import mu_r3bp


# === Parameters === #
mu = mu_r3bp
J2 = J2_neptune
R_nd = R_neptune_meters / a_triton_meters
J2_eff = J2 * R_nd**2
omega2 = 1 + 3 * J2_eff

# === Effective potential derivatives === #
def compute_second_derivatives(x, y, mu, J2_eff, h=1e-6):
    def U_eff(x, y):
        r1 = np.sqrt((x + mu)**2 + y**2)
        r2 = np.sqrt((x - 1 + mu)**2 + y**2)
        return (
            0.5 * omega2 * (x**2 + y**2)
            + (1 - mu) / r1
            + mu / r2
            + (1 - mu) * J2_eff / r1**3
        )

    Uxx = (U_eff(x + h, y) - 2 * U_eff(x, y) + U_eff(x - h, y)) / h**2
    Uyy = (U_eff(x, y + h) - 2 * U_eff(x, y) + U_eff(x, y - h)) / h**2
    Uxy = (
        U_eff(x + h, y + h) - U_eff(x + h, y - h)
        - U_eff(x - h, y + h) + U_eff(x - h, y - h)
    ) / (4 * h**2)

    return Uxx, Uxy, Uyy

# === Linearized matrix and eigenvalue analysis === #
def evaluate_stability(x, y, mu, J2_eff):
    Uxx, Uxy, Uyy = compute_second_derivatives(x, y, mu, J2_eff)

    A = np.array([
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [Uxx, Uxy, 0, 2],
        [Uxy, Uyy, -2, 0],
    ])

    eigenvalues = eigvals(A)
    stability = "Stable" if np.all(np.real(eigenvalues) <= 0) else "Unstable"
    return eigenvalues, stability

# === Main Execution === #
if __name__ == "__main__":
    points = get_lagrange_points_and_constants()
    print("\n📊 Stability Analysis with Neptune’s J2 Perturbation:\n")

    for label in ["L1", "L2", "L3", "L4", "L5"]:
        x, y = points[label]
        eigs, status = evaluate_stability(x, y, mu, J2_eff)
        print(f"🔹 {label} at ({x:.6f}, {y:.6f}): {status}")
        for i, val in enumerate(eigs):
            print(f"    λ{i+1} = {val.real:.6e} {'+' if val.imag >= 0 else '-'} {abs(val.imag):.6e}j")
        print()
