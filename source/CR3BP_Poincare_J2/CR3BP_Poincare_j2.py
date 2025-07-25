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

# === Path Setup ===
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'neptune_triton_models', 'r3bp_with_perturbations')))
from constants import M_neptune, M_triton, R_neptune_meters, J2_neptune, a_triton_meters
from r3bp_calculations_02 import gradient_U

# === System Parameters ===
mu = M_triton / (M_neptune + M_triton)  # Dimensionless
R_nd = R_neptune_meters / a_triton_meters
J2 = J2_neptune

# === User-Defined Parameters ===
params = {
    'C': 2.999,
    'mu': mu,
    'J2': J2,
    'R_nd': R_nd,
    'J2_enabled': False,  # Toggle perturbations
    'integrator': 'DOP853',  # or 'LSODA'
    'tlim': 10000.0,
    'dt': 0.01,
    'folder': "poincare_output",
    'plot': True,
    'combined_plot': True
}

os.makedirs(params['folder'], exist_ok=True)

# === Dynamics ===
def gradient_U_oblate(x, y, z, mu, J2, R):
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x - (1 - mu))**2 + y**2 + z**2)
    A = (3 * J2 * R**2) / 2
    f = A * (1 - 5 * z**2 / r1**2) / r1**5

    Omega_x = x - (1 - mu)*(x + mu)/r1**3 - mu*(x - (1 - mu))/r2**3 + f * (x + mu)
    Omega_y = y - (1 - mu)*y/r1**3 - mu*y/r2**3 + f * y
    Omega_z = -(1 - mu)*z/r1**3 - mu*z/r2**3 + A * z * (3 - 5 * z**2 / r1**2) / r1**5
    return Omega_x, Omega_y, Omega_z

def equations_cr3bp(t, state, mu):
    x, y, z, vx, vy, vz = state
    Ux, Uy, Uz = gradient_U(x, y, z, mu)
    return [vx, vy, vz, Ux + 2 * vy, Uy - 2 * vx, Uz]

def equations_oblate(t, state, mu, J2, R):
    x, y, z, vx, vy, vz = state
    Ux, Uy, Uz = gradient_U_oblate(x, y, z, mu, J2, R)
    return [vx, vy, vz, Ux + 2 * vy, Uy - 2 * vx, Uz]

# === Events ===
def y_cross_event(t, f):
    return f[1]
y_cross_event.direction = 1
y_cross_event.terminal = False

def escape_or_collision_event(t, f):
    x, y, z = f[0], f[1], f[2]
    R1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    R2 = np.sqrt((x - (1 - mu))**2 + y**2 + z**2)
    R = np.sqrt(x**2 + y**2 + z**2)

    r1_min = 3000 / a_triton_meters
    r2_min = 1000 / a_triton_meters
    r_max = 5.0

    if R1 <= r1_min:
        return -1  # M1 collision
    if R2 <= r2_min:
        return -2  # M2 collision
    if R >= r_max:
        return 1   # Escape
    return 10  # No event

escape_or_collision_event.terminal = True
escape_or_collision_event.direction = 0

# === Poincaré Generator ===
def generate_poincare(x0):
    C = params['C']
    mu = params['mu']
    mu_star = 1 - mu
    y0 = 0.0
    z0 = 0.0

    R1 = np.sqrt((x0 + mu)**2)
    R2 = np.sqrt((x0 - mu_star)**2)
    arg = -C + x0**2 + 2 * (mu_star / R1 + mu / R2)
    if arg < 0:
        return f"Skipping x0 = {x0} (imaginary vy0)"

    vy0 = np.sqrt(arg)
    initial_state = [x0, y0, z0, 0.0, vy0, 0.0]

    if params['J2_enabled']:
        rhs = lambda t, f: equations_oblate(t, f, mu, params['J2'], params['R_nd'])
    else:
        rhs = lambda t, f: equations_cr3bp(t, f, mu)

    sol = solve_ivp(
        rhs,
        [0, params['tlim']],
        initial_state,
        t_eval=np.arange(0, params['tlim'], params['dt']),
        events=[y_cross_event, escape_or_collision_event],
        method=params['integrator'],
        rtol=1e-10,
        atol=1e-12,
        max_step=1.0
    )

    if sol.status == 1 and len(sol.t_events[1]) > 0:
        t_event = sol.t_events[1][0]
        f_event = sol.y_events[1][0]
        reason = escape_or_collision_event(t_event, f_event)
        reason_str = { -1: "Collision with M1", -2: "Collision with M2", 1: "Escape" }.get(reason, "Unknown")
        with open(os.path.join(params['folder'], "collisions_and_escapes.log"), "a") as f:
            f.write(f"x0 = {x0}, C = {C}, reason: {reason_str}, t = {t_event:.3f}\n")
        return f"x0 = {x0} -> Excluded: {reason_str}"

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
        return f"x0 = {x0}, crossings = 0"

    data = np.array(poincare_points)
    times = np.array(times)

    fname = f"PY-C{C}_Xi{x0}.dat"
    np.savetxt(os.path.join(params['folder'], fname), data)

    if params['plot']:
        cmap = cm.get_cmap("viridis")
        norm = Normalize(vmin=0, vmax=params['tlim'])
        colors = cmap(norm(times))
        plt.figure(figsize=(8, 6))
        plt.scatter(data[:, 0], data[:, 3], c=colors, s=1)
        plt.xlabel(r'$x$')
        plt.ylabel(r'$\dot{x}$')
        plt.title(f"Poincar\u00e9 Section\nC = {C}, x0 = {x0}")
        plt.colorbar(label='Time')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(params['folder'], f"Poincare_C{C}_x0{x0}.png"), dpi=300)
        plt.close()

    return f"x0 = {x0}, crossings = {len(data)}"

# === Batch Runner ===
if __name__ == "__main__":
    print(f"[Info] R_nd: {R_nd:.6f}, mu: {mu:.8f}")
    print("Start:", datetime.now())
    start = timeit.default_timer()

    x0_values = np.arange(-1.0, 1.0, 0.1)

    with mp.Pool(mp.cpu_count()) as pool:
        results = pool.map(generate_poincare, x0_values)

    for res in results:
        print(res)

    if params['combined_plot']:
        cmap = cm.get_cmap("plasma")
        norm = Normalize(vmin=min(x0_values), vmax=max(x0_values))
        plt.figure(figsize=(8, 6))
        for x0 in x0_values:
            path = os.path.join(params['folder'], f"PY-C{params['C']}_Xi{x0}.dat")
            if os.path.exists(path):
                data = np.loadtxt(path)
                if data.ndim == 1:
                    data = data[np.newaxis, :]
                plt.scatter(data[:, 0], data[:, 3], color=cmap(norm(x0)), s=1, alpha=0.6)
        plt.xlabel(r'$x$')
        plt.ylabel(r'$\dot{x}$')
        plt.title(f"Combined Poincar\u00e9 Section (C = {params['C']})")
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        plt.colorbar(sm, label=r'$x_0$')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(params['folder'], f"Combined_Poincare_C{params['C']}.png"), dpi=300)
        plt.close()

    stop = timeit.default_timer()
    print("End:", datetime.now())
    print(f"Total runtime: {stop - start:.2f} seconds")






