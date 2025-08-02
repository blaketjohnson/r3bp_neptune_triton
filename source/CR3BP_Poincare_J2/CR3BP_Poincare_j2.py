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
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from datetime import datetime
import os
from constants import *
from parameters import *

# === Derived Parameters ===
mu = M_triton / (M_neptune + M_triton)
R_nd = R_neptune_meters / a_triton_meters
r1_min = R_neptune_meters / a_triton_meters  # Neptune collision radius (ND)
r2_min = R_trition_meters / a_triton_meters  # Triton collision radius (ND)

# === Escape Threshold (ND) ===
hill_radius_nd = (a_triton_meters * (mu / 3)**(1/3)) * 3 / a_triton_meters
soi_radius_nd = SOI_neptune * 1000 / a_triton_meters
r_max = soi_radius_nd if use_SOI_as_escape else hill_radius_nd

# === Non-dimensionalization Utilities ===
def time_to_nondimensional(t_sec):
    return t_sec / T_triton_seconds

# === Gradient of Potential Functions ===
def gradient_U(x, y, z):
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x - (1 - mu))**2 + y**2 + z**2)
    Ox = x - (1 - mu)*(x + mu)/r1**3 - mu*(x - (1 - mu))/r2**3
    Oy = y - (1 - mu)*y/r1**3 - mu*y/r2**3
    Oz = -(1 - mu)*z/r1**3 - mu*z/r2**3
    return Ox, Oy, Oz

def gradient_U_oblate(x, y, z):
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x - (1 - mu))**2 + y**2 + z**2)
    A = (3 * J2_neptune * R_nd**2) / 2
    f = A * (1 - 5 * z**2 / r1**2) / r1**5
    Ox = x - (1 - mu)*(x + mu)/r1**3 - mu*(x - (1 - mu))/r2**3 + f * (x + mu)
    Oy = y - (1 - mu)*y/r1**3 - mu*y/r2**3 + f * y
    Oz = -(1 - mu)*z/r1**3 - mu*z/r2**3 + A * z * (3 - 5 * z**2 / r1**2) / r1**5
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
    x, y, z = f[0], f[1], f[2]
    R1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    R2 = np.sqrt((x - (1 - mu))**2 + y**2 + z**2)
    R_total = np.sqrt(x**2 + y**2 + z**2)
    if R1 <= r1_min:
        return -1  # Neptune collision
    if R2 <= r2_min:
        return -2  # Triton collision
    if R_total >= r_max:
        return 1   # Escape
    return 10      # No event
escape_or_collision_event.terminal = True
escape_or_collision_event.direction = 0

# === Output Setup ===
data_folder = "Poincare_data"
os.makedirs(data_folder, exist_ok=True)
log_path = os.path.join(data_folder, "run_log.txt")
log_file = open(log_path, "w")

def generate_poincare(args):
    C, x0 = args
    y0, z0 = YI, 0.0
    R1 = np.sqrt((x0 + mu)**2 + y0**2)
    R2 = np.sqrt((x0 - (1 - mu))**2 + y0**2)
    arg = -C + x0**2 + y0**2 + 2 * ((1 - mu) / R1 + mu / R2)

    if arg < 0:
        return f"Skipping x0 = {x0:.5f} (imaginary vy0)"
    vy0 = np.sqrt(arg)
    initial_state = [x0, y0, z0, 0.0, vy0, 0.0]

    sol = solve_ivp(
        equations,
        [0, time_to_nondimensional(tlim_sec)],
        initial_state,
        t_eval=np.arange(0, time_to_nondimensional(tlim_sec), time_to_nondimensional(dt_sec)),
        events=[y_cross_event, escape_or_collision_event],
        method='DOP853',
        rtol=1e-10,
        atol=1e-12
    )

    if sol.status == 1 and len(sol.t_events[1]) > 0:
        t_event = sol.t_events[1][0]
        f_event = sol.y_events[1][0]
        reason_code = escape_or_collision_event(t_event, f_event)
        reason = { -1: "Neptune Collision", -2: "Triton Collision", 1: "Escape" }.get(reason_code, "Other")
        with open(os.path.join(data_folder, "collisions_and_escapes.log"), "a") as f:
            f.write(f"x0 = {x0:.5f}, C = {C}, Reason: {reason}, t = {t_event:.4f}\n")
        return f"x0 = {x0:.5f}, C = {C:.5f} -> Excluded ({reason})"

    crossings, times = [], []
    for i in range(1, sol.y.shape[1]):
        if sol.y[1, i-1] * sol.y[1, i] < 0 and sol.y[1, i] > 0:
            xm = 0.5 * (sol.y[:, i] + sol.y[:, i-1])
            tm = 0.5 * (sol.t[i] + sol.t[i-1])
            crossings.append(xm)
            times.append(tm)

    if not crossings:
        return f"x0 = {x0:.5f}, C = {C:.5f}, crossings = 0"

    data = np.column_stack((np.array(crossings), np.array(times)))
    fname = f"PY-C{C:.5f}_Xi{x0:.5f}.dat"
    np.savetxt(os.path.join(data_folder, fname), data)
    return f"x0 = {x0:.5f}, C = {C:.5f}, crossings = {len(data)}"

if __name__ == "__main__":
    print(f"[Info] mu = {mu:.8f}, R_nd = {R_nd:.6f}, Hill Radius ND = {r_max:.4f}")
    print("Start:", datetime.now())
    start = timeit.default_timer()

    x0_values = np.arange(XI, XF + DX/2, DX)
    print(f"Initial y0 = {YI}")
    Cs = np.arange(C0, CF + dC/2, dC)
    input_pairs = [(C, x0) for C in Cs for x0 in x0_values]

    with mp.Pool(mp.cpu_count()) as pool:
        results = pool.map(generate_poincare, input_pairs)

    for res in results:
        print(res)
        log_file.write(res + "\n")

    stop = timeit.default_timer()
    print("End:", datetime.now())
    print(f"Total runtime: {stop - start:.2f} seconds")
    log_file.write(f"Total runtime: {stop - start:.2f} seconds\n")
    log_file.close()






