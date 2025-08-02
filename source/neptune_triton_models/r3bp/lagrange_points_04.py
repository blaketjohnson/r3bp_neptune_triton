"""
lagrange_points_02.py

Computes Lagrange point locations and their Jacobi constants
for the Neptune–Triton system in the standard CR3BP (no J2 perturbations).

References:
- Murray & Dermott (1999), Solar System Dynamics, Ch. 3
- Gradient of the effective potential: ∇U = 0 (Murray p. 67, Eq. 3.23)

Usage:
    python lagrange_points_02.py
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar, root
from r3bp_calculations_02 import mu_r3bp as mu

# === Primary positions in rotating frame === #
neptune = (-mu, 0)
triton = (1 - mu, 0)

# === Potential derivative along x-axis === #
def dU_dx(x, mu):
    r1 = np.abs(x + mu)
    r2 = np.abs(x - 1 + mu)
    return (
        x
        - (1 - mu) * (x + mu) / r1**3
        - mu * (x - 1 + mu) / r2**3
    )

# === Full gradient of the effective potential === #
def grad_U(xy, mu):
    x, y = xy
    r1 = np.sqrt((x + mu)**2 + y**2)
    r2 = np.sqrt((x - 1 + mu)**2 + y**2)
    dUx = (
        x
        - (1 - mu) * (x + mu) / r1**3
        - mu * (x - 1 + mu) / r2**3
    )
    dUy = (
        y
        - (1 - mu) * y / r1**3
        - mu * y / r2**3
    )
    return np.array([dUx, dUy])

# === Jacobi constant === #
def jacobi_constant(x, y, mu):
    r1 = np.sqrt((x + mu)**2 + y**2)
    r2 = np.sqrt((x - 1 + mu)**2 + y**2)
    U = (
        0.5 * (x**2 + y**2)
        + (1 - mu) / r1
        + mu / r2
    )
    return 2 * U

# === Lagrange point finders === #
def find_L_point(bracket, mu, label="L-point"):
    a, b = bracket
    fa, fb = dU_dx(a, mu), dU_dx(b, mu)
    if fa * fb > 0:
        raise ValueError(f"⚠️ {label} bracket does not straddle a root.")
    sol = root_scalar(dU_dx, args=(mu,), bracket=(a, b), method='brentq')
    if not sol.converged:
        raise RuntimeError(f"❌ {label} root not found.")
    print(f"✅ {label} found at x = {sol.root:.6f}")
    return sol.root

def find_L4_L5(mu):
    sqrt3over2 = np.sqrt(3) / 2
    x0 = 0.5 - mu
    res4 = root(grad_U, x0=[x0, +sqrt3over2], args=(mu,), method='hybr')
    res5 = root(grad_U, x0=[x0, -sqrt3over2], args=(mu,), method='hybr')
    if not res4.success or not res5.success:
        raise RuntimeError("❌ L4 or L5 root-finding failed.")
    print(f"✅ L4 found at ({res4.x[0]:.6f}, {res4.x[1]:.6f})")
    print(f"✅ L5 found at ({res5.x[0]:.6f}, {res5.x[1]:.6f})")
    return res4.x, res5.x

# === Unified calculation === #
def get_lagrange_points_and_constants():
    x_L1 = find_L_point((0.3, 0.98), mu, "L1")
    x_L2 = find_L_point((1.0, 1.2), mu, "L2")
    x_L3 = find_L_point((-1.5, -0.8), mu, "L3")
    L4, L5 = find_L4_L5(mu)

    L1 = (x_L1, 0)
    L2 = (x_L2, 0)
    L3 = (x_L3, 0)
    L4 = tuple(L4)
    L5 = tuple(L5)

    return {
        "L1": L1, "L2": L2, "L3": L3, "L4": L4, "L5": L5,
        "C_L1": jacobi_constant(*L1, mu),
        "C_L2": jacobi_constant(*L2, mu),
        "C_L3": jacobi_constant(*L3, mu),
        "C_L4": jacobi_constant(*L4, mu),
        "C_L5": jacobi_constant(*L5, mu)
    }

# === Run if main === #
if __name__ == "__main__":
    points = get_lagrange_points_and_constants()

    print("\n📐 Jacobi Constants (CR3BP, no J2):")
    for label in ["C_L1", "C_L2", "C_L3", "C_L4", "C_L5"]:
        print(f"{label} = {points[label]:.6f}")

    # Plot results
    plt.figure(figsize=(8, 6))
    plt.plot(*neptune, 'o', label='Neptune', color='blue')
    plt.plot(*triton, 'o', label='Triton', color='orange')

    for label in ["L1", "L2", "L3", "L4", "L5"]:
        x, y = points[label]
        plt.plot(x, y, 'k*', markersize=10)
        plt.text(x + 0.02, y + 0.02, label, fontsize=9)

    plt.title("Lagrange Points - Neptune–Triton (CR3BP)")
    plt.xlabel("x (nondimensional)")
    plt.ylabel("y (nondimensional)")
    plt.grid(True)
    plt.axis('equal')
    plt.legend()
    plt.tight_layout()
    plt.show()
