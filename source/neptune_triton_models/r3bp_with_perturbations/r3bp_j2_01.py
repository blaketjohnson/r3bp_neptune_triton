"""
r3bp_j2_01.py

Models the CR3BP for Neptune–Triton with Neptune's oblateness (J2) included.

References:
- Murray & Dermott, Solar System Dynamics
  - J2 Perturbation: Eq. (4.104), p. 151
  - CR3BP EOM: Eq. (3.19–3.21), p. 68
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import sys

# === Imports === #
from constants import J2_neptune, R_neptune_meters
from non_dimensionalizer import mu_nd, to_nondim_length

# === Nondimensional Setup === #
mu = mu_nd
R_nd = to_nondim_length(R_neptune_meters)

print(f"[Info] R_nd (Neptune radius / Triton SMA): {R_nd:.6f}")
print(f"[Info] mu_nd (mass ratio): {mu:.8f}")

# === Effective Potential and Gradient === #
def U_oblate(x, y, z, mu, J2, R):
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x + mu - 1)**2 + y**2 + z**2)
    U_classical = 0.5 * (x**2 + y**2 + z**2) + (1 - mu) / r1 + mu / r2
    U_J2 = - (1 - mu) * J2 * R**2 / (2 * r1**3) * (1 - 3 * z**2 / r1**2)
    return U_classical + U_J2

def gradient_U_oblate(x, y, z, mu, J2, R, h=1e-6):
    Ux = (U_oblate(x + h, y, z, mu, J2, R) - U_oblate(x - h, y, z, mu, J2, R)) / (2 * h)
    Uy = (U_oblate(x, y + h, z, mu, J2, R) - U_oblate(x, y - h, z, mu, J2, R)) / (2 * h)
    Uz = (U_oblate(x, y, z + h, mu, J2, R) - U_oblate(x, y, z - h, mu, J2, R)) / (2 * h)
    return np.array([Ux, Uy, Uz])

def equations_oblate(t, state, mu, J2, R):
    x, y, z, vx, vy, vz = state
    Ux, Uy, Uz = gradient_U_oblate(x, y, z, mu, J2, R)
    ax = Ux + 2 * vy
    ay = Uy - 2 * vx
    az = Uz
    return [vx, vy, vz, ax, ay, az]

# === Simulation Routines === #
def run_oblateness_simulation():
    state0 = [0.5, 0.0, 0.0, 0.0, 0.1, 0.0]
    t_span = (0, 100)
    t_eval = np.linspace(*t_span, 1000)

    sol = solve_ivp(
        equations_oblate,
        t_span,
        state0,
        t_eval=t_eval,
        args=(mu, J2_neptune, R_nd),
        method='RK45'
    )
    return sol

def plot_oblate_orbit(sol):
    x, y = sol.y[0], sol.y[1]

    plt.figure(figsize=(8, 6))
    plt.plot(x, y, label='Trajectory with J2', color='teal')
    plt.scatter([-mu], [0], label='Neptune', color='blue')
    plt.scatter([1 - mu], [0], label='Triton', color='brown')
    plt.title("CR3BP with Neptune Oblateness (J2)")
    plt.xlabel("x (nondimensional)")
    plt.ylabel("y (nondimensional)")
    plt.axis('equal')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def run_oblateness_simulation_3d():
    state0 = [0.5, 0.0, 0.01, 0.0, 1.0, 0.001]
    t_span = (0, 100)
    t_eval = np.linspace(*t_span, 2000)

    sol = solve_ivp(
        equations_oblate,
        t_span,
        state0,
        t_eval=t_eval,
        args=(mu, J2_neptune, R_nd),
        method='RK45'
    )
    return sol

def plot_oblate_orbit_3d(sol):
    from mpl_toolkits.mplot3d import Axes3D

    x, y, z = sol.y[0], sol.y[1], sol.y[2]

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    ax.plot(x, y, z, label='Trajectory with J2 (3D)', color='teal')
    ax.scatter([-mu], [0], [0], label='Neptune', color='blue')
    ax.scatter([1 - mu], [0], [0], label='Triton', color='brown')

    ax.set_title("CR3BP with Neptune Oblateness (J2) in 3D")
    ax.set_xlabel("x (nondimensional)")
    ax.set_ylabel("y (nondimensional)")
    ax.set_zlabel("z (nondimensional)")
    ax.legend()
    plt.tight_layout()
    plt.show()

# === Main Execution === #
if __name__ == "__main__":
    sol2d = run_oblateness_simulation()
    plot_oblate_orbit(sol2d)

    sol3d = run_oblateness_simulation_3d()
    plot_oblate_orbit_3d(sol3d)


    

