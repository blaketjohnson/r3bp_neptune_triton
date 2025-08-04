"""
CR3BP_Poincare_J2.py
---------------------------------------
Generates Poincaré surface-of-section data files for the classical
Circular Restricted Three-Body Problem (CR3BP) in the Neptune–Triton system,
with optional inclusion of Neptune's J₂ oblateness perturbation.

Features:
    • Supports two mapping modes:
        - "global"  : Fixed-step integration with midpoint interpolation
        - "highres" : Event-detection for exact y=0 upward crossings
    • Optional inclusion of Neptune's J₂ perturbation via parameters.py
    • Runtime termination for trajectories that escape (Hill radius) or collide
      with Neptune/Triton (user-defined minimum approach distance)
    • Separate output folders for each mapping mode
    • Mode-specific .dat filenames for clarity and to prevent overwrites

Output:
    • .dat files containing Poincaré section crossing data
    • Files are saved in:
        Poincare_data_global/  or  Poincare_data_highres/
      depending on mapping_mode
    • Filenames include mode tag, e.g.:
        PY-C3.00000Xi0.10000_global.dat

Usage:
    1. Configure run settings in `parameters.py`:
        - Jacobi constant range
        - x₀ sweep range
        - Mapping mode ("global" or "highres")
        - Minimum safe distances from Neptune/Triton
        - J₂ perturbation toggle
    2. Run the solver:
        $ python CR3BP_Poincare_classical.py
    3. Resulting .dat files can be plotted with `plot_CR3BP_mapper.py`

Author:
    Blake T. Johnson
    Thesis Project
    (c) 2025
"""


import os
import timeit
import numpy as np
from datetime import datetime
from scipy.integrate import solve_ivp
from constants import *       # Physical & orbital constants
from parameters import *      # User-specified run settings

# --- Collision radii (ND) ---
r1_min_km = R_neptune + min_distance_neptune_km
r2_min_km = R_triton + min_distance_triton_km
r1_min = r1_min_km / a_triton_km
r2_min = r2_min_km / a_triton_km

R_nd = R_neptune / a_triton_km

# --- Derived constants ---
mu = M_triton / (M_neptune + M_triton)
mu_star = 1.0 - mu

# --- Hill Radius Calculation (escape distance) ---
def neptune_hill_radius():
    """Compute Neptune's Hill radius relative to the Sun in km and ND."""
    R_H_km = a_neptune_km * (M_neptune / (3 * M_sun)) ** (1/3)
    R_H_ND = R_H_km / a_triton_km
    return R_H_km, R_H_ND

R_H_km, R_H_ND = neptune_hill_radius()
r_max = R_H_ND  # Escape radius in ND

# --- Accelerations ---
def accel_classical(x, y, z, r1, r2):
    '''Equations of motion for the classical CR3BP without J2 perturbation.'''
    ax = x - (mu_star*(x+mu))/r1**3 - (mu*(x - mu_star))/r2**3
    ay = y - (mu_star*y)/r1**3 - (mu*y)/r2**3
    az = -(mu_star*z)/r1**3 - (mu*z)/r2**3
    return ax, ay, az

def accel_j2(x, y, z, r1, r2):
    '''Equations of motion for the CR3BP with J2 perturbation.'''
    A = (3 * J2_neptune * R_nd**2) / 2
    f = A * (1 - 5 * z**2 / r1**2) / r1**5
    ax = x - (mu_star*(x+mu))/r1**3 - (mu*(x - mu_star))/r2**3 + f * (x+mu)
    ay = y - (mu_star*y)/r1**3 - (mu*y)/r2**3 + f * y
    az = -(mu_star*z)/r1**3 - (mu*z)/r2**3 + A * z * (3 - 5 * z**2 / r1**2) / r1**5
    return ax, ay, az

def cr3bp_equations(t, state):
    """Runs the appropriate Equation of motion based on the J2_enabled flag."""
    x, y, z, vx, vy, vz = state
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x - mu_star)**2 + y**2 + z**2)
    if J2_enabled:
        ax, ay, az = accel_j2(x, y, z, r1, r2)
    else:
        ax, ay, az = accel_classical(x, y, z, r1, r2)
    return [vx, vy, vz, ax + 2*vy, ay - 2*vx, az]


# --- Initial condition collision check ---
def check_initial_conditions(x0, y0, z0):
    R1_0 = np.sqrt((x0 + mu)**2 + y0**2 + z0**2)
    R2_0 = np.sqrt((x0 - mu_star)**2 + y0**2 + z0**2)
    if R1_0 <= r1_min or R2_0 <= r2_min:
        return False
    return True

# --- Runtime escape/collision event ---
def escape_or_collision_event(t, state):
    """Stop integration if escape or collision occurs."""
    x, y, z = state[0], state[1], state[2]
    R1 = np.sqrt((x + mu)**2 + y**2 + z**2)       # Neptune distance
    R2 = np.sqrt((x - mu_star)**2 + y**2 + z**2)  # Triton distance
    R_total = np.sqrt(x**2 + y**2 + z**2)         # Distance from barycenter
    # Check conditions
    if R1 <= r1_min:    # Neptune collision
        return 0
    if R2 <= r2_min:    # Triton collision
        return 0
    if R_total >= r_max:  # Escape
        return 0

    return 1  # No event yet

# Configure event properties
escape_or_collision_event.terminal = True
escape_or_collision_event.direction = -1

# --- Poincaré section generation ---
def y_cross_event(t, state):
    return state[1]  # y coordinate
y_cross_event.direction = 1
y_cross_event.terminal = False

# --- Main integration loop for one trajectory ---
def generate_poincare(C, x0):
    # Choose initial y0 based on mapping mode
    if mapping_mode.lower() == "global":
        y0 = 0.0
    else:
        y0 = YI
    z0 = 0.0


    # Skip if collision at start
    if not check_initial_conditions(x0, y0, z0):
        return f"Skip x0={x0:.5f} - initial collision"

    # vy0 from Jacobi constant
    R1_0 = np.sqrt((x0 + mu)**2 + y0**2 + z0**2)
    R2_0 = np.sqrt((x0 - mu_star)**2 + y0**2 + z0**2)
    arg = -C + x0**2 + y0**2 + 2*((mu_star)/R1_0 + mu/R2_0)
    if arg < 0:
        return f"Skip x0={x0:.5f} - imaginary vy0"
    vy0 = np.sqrt(arg)

    # Initial state
    state0 = [x0, y0, z0, 0.0, vy0, 0.0]
    fname = f"PY-C{C:.5f}Xi{round(x0, 5)}_{mapping_mode.lower()}.dat"
    fpath = os.path.join(output_folder, fname)
    if os.path.exists(fpath): os.remove(fpath)

    # --- HIGHRES mode ---
    if mapping_mode.lower() == "highres":
        sol = solve_ivp(
            cr3bp_equations,
            [0, tlim_sec],
            state0,
            events=[y_cross_event, escape_or_collision_event],
            method='DOP853',
            rtol=1e-12,
            atol=1e-12
        )
        for crossing in sol.y_events[0]:
            line = list(crossing) + [x0, vy0, C]
            with open(fpath, "a") as f:
                np.savetxt(f, [line], fmt="%.14e")

    # --- GLOBAL mode ---
    else:
        t_eval = np.arange(0, tlim_sec + dt_sec, dt_sec)
        sol = solve_ivp(
            cr3bp_equations,
            [0, tlim_sec],
            state0,
            t_eval=t_eval,
            events=[escape_or_collision_event],
            method='DOP853',
            rtol=1e-12,
            atol=1e-12
        )
        y_vals = sol.y[1]
        for i in range(1, len(y_vals)):
            if y_vals[i-1] <= 0 and y_vals[i] > 0:  # upward crossing
                crossing_state = 0.5 * (sol.y[:, i-1] + sol.y[:, i])
                line = list(crossing_state) + [x0, vy0, C]
                with open(fpath, "a") as f:
                    np.savetxt(f, [line], fmt="%.14e")

    return f"x0={x0:.5f} - file saved"

# --- Run all ---
if __name__ == "__main__":
    print("CR3BP Poincaré Map Generator")
    print(f"Escape radius (Hill): {R_H_km:.2f} km ({R_H_ND:.5f} ND)")
    print(f"Neptune collision limit: {r1_min_km:.1f} km ({r1_min:.5f} ND)")
    print(f"Triton  collision limit: {r2_min_km:.1f} km ({r2_min:.5f} ND)")
    if mapping_mode.lower() == "global":
        print("Mapping mode: GLOBAL (fixed-step, midpoint interpolation, y0=0.0)")
    elif mapping_mode.lower() == "highres":
        print(f"Mapping mode: HIGHRES (event detection, y0={YI})")
    else:
        print(f"Mapping mode: {mapping_mode.upper()} (unrecognized, defaulting to global)")
    print("Start:", datetime.now())

    output_folder = f"Poincare_data_{mapping_mode.lower()}"
    os.makedirs(output_folder, exist_ok=True)


    start_time = timeit.default_timer()
    Cs = np.arange(C0, CF + dC/2, dC)
    x0_values = np.arange(XI, XF + DX/2, DX)
    for C in Cs:
        for x0 in x0_values:
            msg = generate_poincare(C, x0)
            print(msg)
    elapsed = timeit.default_timer() - start_time
    print(f"Done in {elapsed/60:.2f} min")
    print("End:", datetime.now())

















