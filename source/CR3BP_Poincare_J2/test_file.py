"""
CR3BP_Poincare_classical.py

Generates Poincaré section .dat files for the classical (non-J2) CR3BP.
Output format matches professor's: x y z xdot ydot zdot x0 vy0 C.

Author: Blake T. Johnson (modernized from Diogo Merguizo Sanchez's style)
"""

import os
import timeit
import numpy as np
from datetime import datetime
from scipy.integrate import solve_ivp
from constants import *
from parameters import *

# --- Derived constants ---
mu = M_triton / (M_neptune + M_triton)
mu_star = 1.0 - mu
r1_min = R_neptune_meters / a_triton_meters  # ND
r2_min = R_trition_meters / a_triton_meters  # ND
r_max = 1.17e8 / a_triton_km  # ND from parameters.py

# --- Event functions ---
def y_cross_event(t, state):
    """Event: crossing y=0 going upward (Poincaré section)."""
    return state[1]
y_cross_event.direction = 1
y_cross_event.terminal = False

def escape_or_collision_event(t, state):
    """Stop integration if collision with primary or escape radius exceeded."""
    x, y, z = state[0], state[1], state[2]
    R1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    R2 = np.sqrt((x - (1 - mu))**2 + y**2 + z**2)
    R_total = np.sqrt(x**2 + y**2 + z**2)
    if R1 <= r1_min or R2 <= r2_min or R_total >= r_max:
        return 0
    return 1
escape_or_collision_event.terminal = True
escape_or_collision_event.direction = -1

# --- Classical CR3BP equations of motion ---
def cr3bp_equations(t, state):
    x, y, z, vx, vy, vz = state
    R1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    R2 = np.sqrt((x - (1 - mu))**2 + y**2 + z**2)
    ax = 2*vy + x - (mu_star*(x + mu))/R1**3 - (mu*(x - mu_star))/R2**3
    ay = -2*vx + y - (mu_star*y)/R1**3 - (mu*y)/R2**3
    az = -(mu_star*z)/R1**3 - (mu*z)/R2**3
    return [vx, vy, vz, ax, ay, az]

# --- Main integration loop ---
def generate_poincare(C, x0):
    """Integrate one trajectory and save Poincaré crossings."""
    y0, z0 = 0.0, 0.0

    # Collision/escape check
    R1_0 = np.sqrt((x0 + mu)**2 + y0**2 + z0**2)
    R2_0 = np.sqrt((x0 - mu_star)**2 + y0**2 + z0**2)
    if R1_0 <= r1_min or R2_0 <= r2_min:
        return f"Skip x0={x0:.5f} - initial collision"

    # vy0 from Jacobi constant
    arg = -C + x0**2 + y0**2 + 2*((mu_star)/R1_0 + mu/R2_0)
    if arg < 0:
        return f"Skip x0={x0:.5f} - imaginary vy0"
    vy0 = np.sqrt(arg)

    # Initial state
    state0 = [x0, y0, z0, 0.0, vy0, 0.0]

    # Integrate
    sol = solve_ivp(
        cr3bp_equations,
        [0, tlim_sec],
        state0,
        events=[y_cross_event, escape_or_collision_event],
        method='DOP853',
        max_step=dt_sec,
        first_step=dt_sec/100.0,
        rtol=1e-12,
        atol=1e-12
    )

    # If only escape/collision, skip
    if sol.status == 1 and len(sol.t_events[1]) > 0:
        return f"Skip x0={x0:.5f} - escape/collision"

    # Save crossings
    fname = f"PY-C{C:.5f}Xi{round(x0, 5)}.dat"
    fpath = os.path.join(output_folder, fname)
    if os.path.exists(fpath):
        os.remove(fpath)

    for k in range(sol.t_events[0].size):
        crossing_state = sol.y_events[0][k]
        line = list(crossing_state) + [x0, vy0, C]
        with open(fpath, "a") as f:
            np.savetxt(f, [line], fmt="%.14e")

    return f"x0={x0:.5f} - {sol.t_events[0].size} crossings"

# --- Run all ---
if __name__ == "__main__":
    print("Classical CR3BP Poincaré Map Generator")
    print("Output format matches professor's .dat files")
    print("Start:", datetime.now())

    # Create output folder
    output_folder = "Poincare_data_CR3BP"
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















