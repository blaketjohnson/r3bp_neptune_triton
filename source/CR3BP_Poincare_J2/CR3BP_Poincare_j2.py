import timeit
import numpy as np
import multiprocessing as mp
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from datetime import datetime
import os

import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'neptune_triton_models', 'r3bp_with_perturbations')))

from functions import print_status
from constants import M_neptune, M_triton, R_neptune_meters, J2_neptune, a_triton_meters

# === System parameters ===
mu = M_triton / (M_neptune + M_triton) # Dimensionless mass parameter
R_nd = R_neptune_meters / a_triton_meters # Neptune radius in dimensionless units
J2 = J2_neptune # J2 coefficient for Neptune

# === Simulation Parameters ===
params = {
    'C': 2.999,
    'mu': mu,
    'J2': J2,
    'R_nd': R_nd,
    'tlim': 5000.0,
    'dt': 0.01,
    'folder': "poincare_output",
    'plot': True,
    'combined_plot': True
}

os.makedirs(params['folder'], exist_ok=True)

# === Worker Function ===
def gradient_U_oblate(x, y, z, mu, J2, R):
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x - (1 - mu))**2 + y**2 + z**2)
    A = (3 * J2 * R**2) / 2
    f = A * (1 - 5 * z**2 / r1**2) / r1**5

    Omega_x = x - (1 - mu)*(x + mu)/r1**3 - mu*(x - (1 - mu))/r2**3 + f * (x + mu)
    Omega_y = y - (1 - mu)*y/r1**3 - mu*y/r2**3 + f * y
    Omega_z = -(1 - mu)*z/r1**3 - mu*z/r2**3 + A * z * (3 - 5 * z**2 / r1**2) / r1**5

    return Omega_x, Omega_y, Omega_z

def equations_oblate(t, state, mu, J2, R):
    x, y, z, vx, vy, vz = state
    Ux, Uy, Uz = gradient_U_oblate(x, y, z, mu, J2, R)
    ax = Ux + 2 * vy
    ay = Uy - 2 * vx
    az = Uz
    return [vx, vy, vz, ax, ay, az]

def generate_poincare(x0):
    C = params['C']
    mu = params['mu']
    mu_star = 1 - mu

    y0 = np.sqrt(3) / 2
    #y0 = 0
    R1 = np.sqrt((x0 + mu)**2 + y0**2)
    R2 = np.sqrt((x0 - mu_star)**2 + y0**2)
    arg = -C + x0**2 + y0**2 + 2.0 * (mu_star / R1 + mu / R2)
    if arg < 0:
        return f"Skipping x0 = {x0:.3f} (no real vy)"

    vy0 = np.sqrt(arg)
    x = np.zeros(6)
    x[0] = x0
    x[1] = y0
    x[4] = vy0

    def y_cross_event(t, f):
        return f[1] - np.sqrt(3)/2  # trigger when y ≈ √3/2
    y_cross_event.direction = 1
    y_cross_event.terminal = False

    sol = solve_ivp(
        lambda t, f: equations_oblate(t, f, mu, params['J2'], params['R_nd']),
        [0, params['tlim']],
        x,
        t_eval=np.arange(0, params['tlim'], params['dt']),
        events=[y_cross_event],
        rtol=1e-7,
        atol=1e-10,
        max_step=1.0
    )

    poincare_points = []
    times = []
    state = sol.y
    for i in range(1, state.shape[1]):
        if state[1, i-1] * state[1, i] < 0 and state[1, i] > 0:
            xm = 0.5 * (state[:, i] + state[:, i-1])
            tm = 0.5 * (sol.t[i] + sol.t[i-1])
            
            # Safeguard against runaway orbits
            if abs(xm[0]) > 3 or abs(xm[3]) > 2:
                continue  # Skip this crossing if x or vx is too extreme
            
            poincare_points.append(xm)
            times.append(tm)


    if not poincare_points:
        return f"x0 = {x0:.3f}, crossings = 0"

    data = np.array(poincare_points)
    times = np.array(times)
    dat_path = os.path.join(params['folder'], f"PY-C{C:.3f}Xi{x0:.3f}.dat")
    np.savetxt(dat_path, data)

    if params['plot']:
        cmap = cm.get_cmap("viridis")
        norm = Normalize(vmin=0, vmax=params['tlim'])
        colors = cmap(norm(times))

        plt.figure(figsize=(8, 6))
        plt.scatter(data[:, 0], data[:, 3], c=colors, s=1)
        plt.xlabel(r'$x$', fontsize=14)
        plt.ylabel(r'$\dot{x}$', fontsize=14)
        plt.title(f"Poincar\'e Section\nC={C:.3f}, x0={x0:.3f}")
        plt.colorbar(label='Time')
        plt.grid(True)
        plt.tight_layout()
        plot_path = os.path.join(params['folder'], f"Poincare_C{C:.3f}_x0{x0:.3f}.png")
        plt.savefig(plot_path, dpi=300)
        plt.close()

    return f"x0 = {x0:.3f}, crossings = {len(data)}"

# === Parallel Execution ===
if __name__ == "__main__":
    print(f"[Info] R_nd: {R_nd:.6f}, mu: {mu:.8f}")
    print("Init:", datetime.now())
    start = timeit.default_timer()

    x0_values = np.round(np.arange(0.490, 0.515, 0.001), 3)


    with mp.Pool(mp.cpu_count()) as pool:
        results = pool.map(generate_poincare, x0_values)

    for res in results:
        print(res)

    if params.get('combined_plot', False):
        cmap_combined = plt.get_cmap("plasma")
        norm_combined = Normalize(vmin=min(x0_values), vmax=max(x0_values))
        plt.figure(figsize=(8, 6))
        for x0 in x0_values:
            dat_path = os.path.join(params['folder'], f"PY-C{params['C']:.3f}Xi{x0:.3f}.dat")
            if os.path.exists(dat_path):
                data = np.loadtxt(dat_path)
                if data.ndim == 1:
                    data = data[np.newaxis, :]
                color = cmap_combined(norm_combined(x0))
                plt.scatter(data[:, 0], data[:, 3], color=color, s=1, alpha=0.6)
        plt.xlabel(r'$x$', fontsize=14)
        plt.ylabel(r'$\dot{x}$', fontsize=14)
        plt.title(f"Combined Poincar\'e Section (C = {params['C']:.2f})")
        sm = plt.cm.ScalarMappable(cmap=cmap_combined, norm=norm_combined)
        sm.set_array([])
        ax = plt.gca()
        plt.colorbar(sm, ax=ax, label=r'$x_0$')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(params['folder'], f"Combined_Poincare_C{params['C']:.2f}.png"), dpi=300)
        plt.close()

    stop = timeit.default_timer()
    print("\nEnd:", datetime.now())
    print(f"Total runtime: {stop - start:.2f} seconds")




