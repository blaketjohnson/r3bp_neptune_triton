"""
test.py

Quick skip-logic test for CR3BP_Poincare_J2 to confirm initial
collision/escape and imaginary vy0 conditions match professor's behavior.

Author: Blake T. Johnson
"""

import numpy as np
from constants import *
from parameters import *

# === Derived parameters ===
mu = M_triton / (M_neptune + M_triton)

# ND distances for collision checks
r1_min = R_neptune_meters / a_triton_meters
r2_min = R_trition_meters / a_triton_meters

# Hill radius (km and ND)
def neptune_hill_radius(M_neptune, M_sun, a_neptune_au, AU_km, a_triton_km):
    a_neptune_km = a_neptune_au * AU_km
    R_H_km = a_neptune_km * (M_neptune / (3 * M_sun)) ** (1/3)
    R_H_ND = R_H_km / a_triton_km
    return R_H_km, R_H_ND

R_H_km, R_H_ND = neptune_hill_radius(
    M_neptune=M_neptune,
    M_sun=1.9885e30,          # kg
    a_neptune_au=30.06896348, # AU
    AU_km=1.495978707e8,      # km
    a_triton_km=a_triton_km
)
r_max = R_H_ND

print(f"\n[Info] mu = {mu:.8f}")
print(f"[Info] r1_min (ND) = {r1_min:.6f}")
print(f"[Info] r2_min (ND) = {r2_min:.6f}")
print(f"[Info] r_max  (ND) = {r_max:.6f}")
print(f"[Info] Hill radius (km) = {R_H_km:.3f}")

# === Sweep through x0 values ===
x0_values = np.arange(XI, XF + DX/2, DX)

print("\n=== Skip Logic Test ===")
for x0 in x0_values:
    y0, z0 = YI, 0.0
    R1_0 = np.sqrt((x0 + mu)**2 + y0**2 + z0**2)
    R2_0 = np.sqrt((x0 - (1 - mu))**2 + y0**2 + z0**2)
    Rtot_0 = np.sqrt(x0**2 + y0**2 + z0**2)

    # Default status
    status = "Integrate"

    # Initial collision/escape checks
    if R1_0 <= r1_min:
        status = "Skip: Initial collision w/ Neptune"
    elif R2_0 <= r2_min:
        status = "Skip: Initial collision w/ Triton"
    elif Rtot_0 >= r_max:
        status = "Skip: Initial escape"

    # Imaginary vy0 check
    else:
        arg = -C0 + x0**2 + y0**2 + 2 * ((1 - mu) / R1_0 + mu / R2_0)
        print(f"x0 = {x0:+.3f} -> arg = {arg:.10f}")
        if arg < 0:
            status = "Skip: Imaginary vy0"


    print(f"x0 = {x0:+.3f} -> {status}")







