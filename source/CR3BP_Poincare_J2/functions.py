"""
functions.py

Core physics, event logic, and utility functions for CR3BP and J2-perturbed models in the Neptune-Triton system.
All dynamics are defined in dimensionless form; this file includes a nondimensionalizer based on `a_triton_meters`.

Author: Blake T. Johnson (modularized)
"""

import numpy as np
from constants import *

# === Non-dimensionalization ===
def dimensional_to_nondimensional(x_km):
    return x_km * 1000 / a_triton_meters  # converts km to ND

def time_to_nondimensional(t_sec):
    return t_sec / T_triton_seconds

def velocity_to_nondimensional(v_mps):
    return v_mps / (a_triton_meters / T_triton_seconds)

# === System Parameters (derived from constants) ===
mu = M_triton / (M_neptune + M_triton)
R_nd = R_neptune_meters / a_triton_meters

# === Gradient of Potential: Unperturbed CR3BP ===
def gradient_U(x, y, z, mu):
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x - (1 - mu))**2 + y**2 + z**2)
    
    Omega_x = x - (1 - mu)*(x + mu)/r1**3 - mu*(x - (1 - mu))/r2**3
    Omega_y = y - (1 - mu)*y/r1**3 - mu*y/r2**3
    Omega_z = -(1 - mu)*z/r1**3 - mu*z/r2**3

    return Omega_x, Omega_y, Omega_z

# === Gradient of Potential: Oblateness-perturbed CR3BP ===
def gradient_U_oblate(x, y, z, mu, J2, R):
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x - (1 - mu))**2 + y**2 + z**2)
    A = (3 * J2 * R**2) / 2
    f = A * (1 - 5 * z**2 / r1**2) / r1**5

    Omega_x = x - (1 - mu)*(x + mu)/r1**3 - mu*(x - (1 - mu))/r2**3 + f * (x + mu)
    Omega_y = y - (1 - mu)*y/r1**3 - mu*y/r2**3 + f * y
    Omega_z = -(1 - mu)*z/r1**3 - mu*z/r2**3 + A * z * (3 - 5 * z**2 / r1**2) / r1**5

    return Omega_x, Omega_y, Omega_z

# === Equations of Motion ===
def equations_cr3bp(t, state, mu):
    x, y, z, vx, vy, vz = state
    Ux, Uy, Uz = gradient_U(x, y, z, mu)
    ax = Ux + 2 * vy
    ay = Uy - 2 * vx
    az = Uz
    return [vx, vy, vz, ax, ay, az]

def equations_oblate(t, state, mu, J2, R):
    x, y, z, vx, vy, vz = state
    Ux, Uy, Uz = gradient_U_oblate(x, y, z, mu, J2, R)
    ax = Ux + 2 * vy
    ay = Uy - 2 * vx
    az = Uz
    return [vx, vy, vz, ax, ay, az]

# === Poincare Section at y = 0 ===
def y_cross_event(t, f):
    return f[1]  # y = 0

y_cross_event.direction = 1

y_cross_event.terminal = False

# === Escape or Collision Event ===
def escape_or_collision_event(t, f):
    x, y, z = f[0], f[1], f[2]
    R1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    R2 = np.sqrt((x - (1 - mu))**2 + y**2 + z**2)
    R_total = np.sqrt(x**2 + y**2 + z**2)

    r1_min = R_neptune_meters / a_triton_meters  # ND
    r2_min = R_trition_meters / a_triton_meters  # ND
    r_max = 5.0  # ND, arbitrary escape limit

    if R1 <= r1_min:
        return -1  # Collision with Neptune
    if R2 <= r2_min:
        return -2  # Collision with Triton
    if R_total >= r_max:
        return 1   # Escape

    return 10  # No event

escape_or_collision_event.terminal = True
escape_or_collision_event.direction = 0

# === Optional: Utility Function for Displaying Info ===
def print_status(x0, C, reason, t):
    print(f"[!] Excluded orbit at x0 = {x0:.4f}, C = {C:.5f}, reason = {reason}, t = {t:.2f}")



