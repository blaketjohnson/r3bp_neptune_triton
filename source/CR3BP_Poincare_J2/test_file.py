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
from functions import *


# === User-Defined Parameters (Dimensional Input) ===
C = 2.999
XI = 0.200  # non-dimensional
XF = 0.215  # non-dimensional
DX = 0.005  # non-dimensional
J2_enabled = False  # Toggle J2 perturbation
plot = True
combined_plot = True

# === Time Parameters (Dimensional Input) ===
tlim_sec = 10000.0  # seconds
dt_sec = 0.0001     # seconds

# === Output Folder ===
output_folder = "poincare_output"
os.makedirs(output_folder, exist_ok=True)

# === Derived Grid ===
x0_values_nd = np.arange(XI, XF, DX)

# === Core Function ===
def generate_poincare(x0_nd):
    x0 = x0_nd
    y0 = 0.0
    z0 = 0.0

    mu_star = 1 - mu
    R1 = np.sqrt((x0 + mu)**2)
    R2 = np.sqrt((x0 - mu_star)**2)
    arg = -C + x0**2 + 2 * (mu_star / R1 + mu / R2)
    if arg < 0:
        return f"Skipping x0 = {x0:.5f} (imaginary vy0)"

    vy0 = np.sqrt(arg)
    initial_state = [x0, y0, z0, 0.0, vy0, 0.0]

    rhs = (lambda t, f: equations_oblate(t, f, mu, J2_neptune, R_nd)) if J2_enabled else \
          (lambda t, f: equations_cr3bp(t, f, mu))

    sol = solve_ivp(
        rhs,
        [0, time_to_nondimensional(tlim_sec)],
        initial_state,
        t_eval=np.arange(0, time_to_nondimensional(tlim_sec), time_to_nondimensional(dt_sec)),
        events=[y_cross_event, escape_or_collision_event],
        method='DOP853',
        rtol=1e-10,
        atol=1e-12,
        max_step=1.0
    )

    if sol.status == 1 and len(sol.t_events[1]) > 0:
        t_event = sol.t_events[1][0]
        f_event = sol.y_events[1][0]
        reason_code = escape_or_collision_event(t_event, f_event)
        reason_str = { -1: "Collision with Neptune", -2: "Collision with Triton", 1: "Escape" }.get(reason_code, "Unknown")
        with open(os.path.join(output_folder, "collisions_and_escapes.log"), "a") as f:
            f.write(f"x0 = {x0:.5f}, C = {C}, reason: {reason_str}, t = {t_event:.3f}\n")
        return f"x0 = {x0:.5f} -> Excluded: {reason_str}"

    poincare_points = []
    times = []
    state = sol.y
    for i in range(1, state.shape[1]):
        if state[1, i-1] * state[1, i] < 0 and state[1, i] > 0:
            xm = 0.5 * (state[:, i] + state[:, i-1])
            tm = 0.5 * (sol.t[i] + sol.t[i-1])
            poincare_points.append(xm)
            times.append(tm)

    if not poincare_points:
        return f"x0 = {x0:.5f}, crossings = 0"

    data = np.array(poincare_points)
    times = np.array(times)
    fname = f"PY-C{C}_Xi{x0:.5f}.dat"
    np.savetxt(os.path.join(output_folder, fname), data)

    if plot:
        cmap = plt.get_cmap("viridis")
        norm = Normalize(vmin=0, vmax=time_to_nondimensional(tlim_sec))
        colors = cmap(norm(times))
        plt.figure(figsize=(8, 6))
        plt.scatter(data[:, 0], data[:, 3], c=colors, s=1)
        plt.xlabel(r'$x$')
        plt.ylabel(r'$\dot{x}$')
        plt.title(f"Poincaré Section\nC = {C}, x0 = {x0:.5f}")
        plt.colorbar(label='Time (ND)')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(output_folder, f"Poincare_C{C}_x0{x0:.5f}.png"), dpi=300)
        plt.close()

    return f"x0 = {x0:.5f}, crossings = {len(data)}"

# === Execution Block ===
if __name__ == "__main__":
    print(f"[Info] R_nd: {R_nd:.6f}, mu: {mu:.8f}")
    print("Start:", datetime.now())
    start = timeit.default_timer()

    with mp.Pool(mp.cpu_count()) as pool:
        results = pool.map(generate_poincare, x0_values_nd)

    for res in results:
        print(res)

    if combined_plot:
        cmap = plt.get_cmap("plasma")
        norm = Normalize(vmin=min(x0_values_nd), vmax=max(x0_values_nd))
        fig, ax = plt.subplots(figsize=(8, 6))
        for x0 in x0_values_nd:
            path = os.path.join(output_folder, f"PY-C{C}_Xi{x0:.5f}.dat")
            if os.path.exists(path):
                data = np.loadtxt(path)
                if data.ndim == 1:
                    data = data[np.newaxis, :]
                ax.scatter(data[:, 0], data[:, 3], color=cmap(norm(x0)), s=1, alpha=0.6)
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$\dot{x}$')
        ax.set_title(f"Combined Poincaré Section (C = {C})")
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        plt.colorbar(sm, ax=ax, label=r'$x_0$ [ND]')
        ax.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(output_folder, f"Combined_Poincare_C{C}.png"), dpi=300)
        plt.close()

    stop = timeit.default_timer()
    print("End:", datetime.now())
    print(f"Total runtime: {stop - start:.2f} seconds")



