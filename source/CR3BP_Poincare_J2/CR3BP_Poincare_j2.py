"""
CR3BP_Poincare_j2.py

Main driver script for generating Poincaré sections in the Neptune-Triton CR3BP system,
with optional inclusion of Neptune's J2 perturbation.

Author: Blake T. Johnson
Summer Thesis Project, Phase 1
"""

import timeit
import numpy as np
import multiprocessing as mp
from scipy.integrate import solve_ivp
from datetime import datetime
import os
from constants import *
from parameters import *

# === Derived Parameters ===
mu = M_triton / (M_neptune + M_triton)
R_nd = R_neptune_meters / a_triton_meters
r1_min = R_neptune_meters / a_triton_meters  # Neptune collision radius (ND)
r2_min = R_trition_meters / a_triton_meters  # Triton collision radius (ND)

# === Hill radius function ===
def neptune_hill_radius(M_neptune, M_sun, a_neptune_au, AU_km, a_triton_km):
    """Compute Neptune's Hill radius relative to the Sun, both in km and ND form."""
    a_neptune_km = a_neptune_au * AU_km
    R_H_km = a_neptune_km * (M_neptune / (3 * M_sun)) ** (1/3)
    R_H_ND = R_H_km / a_triton_km
    return R_H_km, R_H_ND

# Compute once at import so workers also have it
R_H_km, R_H_ND = neptune_hill_radius(
    M_neptune=M_neptune,
    M_sun=1.9885e30,          # kg
    a_neptune_au=30.06896348, # AU
    AU_km=1.495978707e8,      # km
    a_triton_km=a_triton_km
)
r_max = R_H_ND  # global ND escape limit

# === Gradient of Potential Functions ===
def gradient_U(x, y, z):
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x - (1 - mu))**2 + y**2 + z**2)
    Ox = x - (1 - mu) * (x + mu) / r1**3 - mu * (x - (1 - mu)) / r2**3
    Oy = y - (1 - mu) * y / r1**3 - mu * y / r2**3
    Oz = -(1 - mu) * z / r1**3 - mu * z / r2**3
    return Ox, Oy, Oz

def gradient_U_oblate(x, y, z):
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x - (1 - mu))**2 + y**2 + z**2)
    A = (3 * J2_neptune * R_nd**2) / 2
    f = A * (1 - 5 * z**2 / r1**2) / r1**5
    Ox = x - (1 - mu) * (x + mu) / r1**3 - mu * (x - (1 - mu)) / r2**3 + f * (x + mu)
    Oy = y - (1 - mu) * y / r1**3 - mu * y / r2**3 + f * y
    Oz = -(1 - mu) * z / r1**3 - mu * z / r2**3 + A * z * (3 - 5 * z**2 / r1**2) / r1**5
    return Ox, Oy, Oz

# === Equations of Motion ===
def equations(t, state):
    x, y, z, vx, vy, vz = state
    if J2_enabled:
        Ux, Uy, Uz = gradient_U_oblate(x, y, z)
    else:
        Ux, Uy, Uz = gradient_U(x, y, z)
    return [vx, vy, vz, Ux + 2 * vy, Uy - 2 * vx, Uz]

# === Event Triggers ===
def y_cross_event(t, f):
    return f[1]
y_cross_event.direction = 1
y_cross_event.terminal = False

def escape_or_collision_event(t, f):
    # Default: continue integration
    status = 1.0  # positive means keep going

    x, y, z = f[0], f[1], f[2]
    R1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    R2 = np.sqrt((x - (1 - mu))**2 + y**2 + z**2)
    R_total = np.sqrt(x**2 + y**2 + z**2)

    # Collision with Neptune
    if R1 <= r1_min:
        status = 0.0

    # Collision with Triton
    elif R2 <= r2_min:
        status = 0.0

    # Escape
    elif R_total >= r_max:
        status = 0.0

    return status  # 1 = continue, 0 = stop

escape_or_collision_event.terminal = True
escape_or_collision_event.direction = -1  # detect decreasing through zero


# === Output Setup ===
data_folder = "Poincare_data"
os.makedirs(data_folder, exist_ok=True)
log_path = os.path.join(data_folder, "run_log.txt")
log_file = open(log_path, "w")

# === Poincaré Map Generation (matches professor's skip logic) ===
def generate_poincare(args):
    C, x0 = args
    y0, z0 = YI, 0.0

    # === Initial ND distances to primaries and barycenter ===
    R1_0 = np.sqrt((x0 + mu)**2 + y0**2 + z0**2)
    R2_0 = np.sqrt((x0 - (1 - mu))**2 + y0**2 + z0**2)
    Rtot_0 = np.sqrt(x0**2 + y0**2 + z0**2)

    # === Skip if initial position is already a collision or escape ===
    if R1_0 <= r1_min or R2_0 <= r2_min or Rtot_0 >= r_max:
        return f"Skipping x0 = {x0:.5f} (initial collision/escape)"

    # === Jacobi constant velocity check (imaginary vy0 case) ===
    R1 = R1_0
    R2 = R2_0
    arg = -C + x0**2 + y0**2 + 2 * ((1 - mu) / R1 + mu / R2)
    if arg < 0:
        return f"Skipping x0 = {x0:.5f} (imaginary vy0)"

    vy0 = np.sqrt(arg)
    initial_state = [x0, y0, z0, 0.0, vy0, 0.0]

    # === Integrate ===
    sol = solve_ivp(
        equations,
        [0, tlim_sec],
        initial_state,
        t_eval=np.arange(0, tlim_sec, dt_sec),
        events=[y_cross_event, escape_or_collision_event],
        method='DOP853',
        first_step=dt_sec / 100.0,
        rtol=1e-10,
        atol=1e-12
    )

    # === Handle escape/collision events ===
    if sol.status == 1 and len(sol.t_events[1]) > 0:
        t_event = sol.t_events[1][0]
        f_event = sol.y_events[1][0]
        reason_code = escape_or_collision_event(t_event, f_event)
        reason = {-1: "Neptune Collision", -2: "Triton Collision", 1: "Escape"}.get(reason_code, "Other")
        with open(os.path.join(data_folder, "collisions_and_escapes.log"), "a") as f:
            f.write(f"x0 = {x0:.5f}, C = {C}, Reason: {reason}, t = {t_event:.4f}\n")
        return f"x0 = {x0:.5f}, C = {C:.5f} -> Excluded ({reason})"

    # === Detect Poincaré crossings ===
    crossings, times = [], []
    for i in range(1, sol.y.shape[1]):
        # Detect y-crossing going upward
        if sol.y[1, i-1] * sol.y[1, i] < 0 and sol.y[1, i] > 0:
            xm = 0.5 * (sol.y[:, i] + sol.y[:, i-1])
            tm = 0.5 * (sol.t[i] + sol.t[i-1])

            # === Velocity sanity filter (ND) ===
            if np.any(np.abs(xm[3:6]) > 10.0):
                continue

            crossings.append(xm)
            times.append(tm)

    # === Skip if no valid crossings remain ===
    if not crossings:
        return f"x0 = {x0:.5f}, C = {C:.5f}, crossings = 0"

    # === Save ND positions and velocities ===
    data = np.column_stack((np.array(crossings), np.array(times)))
    fname = f"PY-C{C:.5f}_Xi{x0:.5f}.dat"
    np.savetxt(os.path.join(data_folder, fname), data)

    return f"x0 = {x0:.5f}, C = {C:.5f}, crossings = {len(data)}"

# === Main Execution ===
if __name__ == "__main__":
    # Print info
    print(f"Neptune Hill radius: {R_H_km:.6f} km")
    print(f"Neptune Hill radius (ND): {R_H_ND:.6f}")
    print(f"[Info] mu = {mu:.8f}, R_nd = {R_nd:.6f}")
    print(f"Initial y0 = {YI}")
    print("\nRunning...")

    start = timeit.default_timer()

    # Prepare work list for multiprocessing
    x0_values = np.arange(XI, XF + DX/2, DX)
    Cs = np.arange(C0, CF + dC/2, dC)
    input_pairs = [(C, x0) for C in Cs for x0 in x0_values]

    # Run in parallel
    with mp.Pool(mp.cpu_count()) as pool:
        results = pool.map(generate_poincare, input_pairs)

    # Output results
    for res in results:
        print(res)
        log_file.write(res + "\n")

    stop = timeit.default_timer()
    runtime = stop - start

    # Runtime formatting
    if runtime < 60:
        runtime_str = f"{runtime:.2f} seconds"
    elif runtime < 3600:
        runtime_str = f"{runtime/60:.2f} minutes"
    else:
        runtime_str = f"{runtime/3600:.2f} hours"

    print("End:", datetime.now())
    print(f"Runtime = {runtime_str}")
    log_file.write(f"Runtime = {runtime_str}\n")
    log_file.close()













