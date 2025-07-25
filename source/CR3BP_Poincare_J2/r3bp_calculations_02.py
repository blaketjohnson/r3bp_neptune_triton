"""
Restricted Three-Body Problem (RTBP) for the Neptune–Triton System

This module calculates the motion of a third body in the RTBP framework.
It uses symbolic gradients of the effective potential to derive and solve the equations of motion.

Author: Blake T. Johnson
Summer Thesis Project, Phase 1
References:
- Curtis, Orbital Mechanics for Engineering Students, Chapter 3
- Murray & Dermott, Solar System Dynamics, p. 66–67
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sympy as sp
from scipy.integrate import solve_ivp
from constants import *

def mass_ratio(M1, M2):
    """
    Calculate the mass ratio for the RTBP.
    Curtis p. 66 Eq. 3.21

    :param M1: Mass of the first body (kg)
    :param M2: Mass of the second body (kg)
    :return: Mass ratio mu
    """
    return M2 / (M1 + M2)

def get_potential_gradient_funcs():
    x, y, z, u = sp.symbols('x y z u')
    r1 = sp.sqrt((x + u)**2 + y**2 + z**2)
    r2 = sp.sqrt((x + u - 1)**2 + y**2 + z**2)
    U_sym = 0.5 * (x**2 + y**2 + z**2) + (1 - u)/r1 + u/r2

    Ux = sp.diff(U_sym, x)
    Uy = sp.diff(U_sym, y)
    Uz = sp.diff(U_sym, z)

    Ux_func = sp.lambdify((x, y, z, u), Ux, 'numpy')
    Uy_func = sp.lambdify((x, y, z, u), Uy, 'numpy')
    Uz_func = sp.lambdify((x, y, z, u), Uz, 'numpy')

    return Ux_func, Uy_func, Uz_func

def gradient_U(x, y, z, mu):
    """
    Compute gradient of effective potential.
    Curtis p. 67 Eq. 3.23–3.25

    :param x, y, z: Coordinates
    :param mu: Mass ratio
    :return: Gradient [Ux, Uy, Uz]
    """
    return np.array([U_x_func(x, y, z, mu),
                     U_y_func(x, y, z, mu),
                     U_z_func(x, y, z, mu)])

def r3bp_equations(t, state, mu):
    """
    RTBP equations of motion in rotating frame.
    Curtis p. 67 Eq. 3.23–3.25

    :param t: Time (s)
    :param state: State vector [x, y, z, vx, vy, vz]
    :param mu: Mass ratio
    :return: Derivatives of state vector
    """
    x, y, z, vx, vy, vz = state
    Ux, Uy, Uz = gradient_U(x, y, z, mu)
    ax = Ux + 2 * vy
    ay = Uy - 2 * vx
    az = Uz
    return [vx, vy, vz, ax, ay, az]

def plot_three_body_orbit(sol, label_orbit="RTBP Orbit", label_primary="Neptune", label_secondary="Triton", label_start="Start", color_orbit='green', color_primary='blue', color_secondary='saddlebrown'):
    """
    Plot the 2D projection of a 3D RTBP orbit.
    Matches style of plot_two_body_orbit.
    
    :param sol: (N x 3) solution array with [x, y, z]
    :param label_orbit: Label for orbit line
    :param label_primary: Label for Neptune
    :param label_secondary: Label for Triton
    :param label_start: Label for initial point
    :param color_orbit: Orbit color
    :param color_primary: Neptune color
    :param color_secondary: Triton color
    """
    plt.figure()
    plt.plot(sol[:, 0], sol[:, 1], color=color_orbit, label=label_orbit)
    plt.scatter(-mu_r3bp, 0, color=color_primary, s=50, label=label_primary)  # Neptune
    plt.scatter(1 - mu_r3bp, 0, color=color_secondary, s=50, label=label_secondary)  # Triton
    plt.scatter(sol[0, 0], sol[0, 1], color=color_orbit, s=30, label=label_start)  # Start
    plt.xlabel("x (non-dimensional)")
    plt.ylabel("y (non-dimensional)")
    plt.legend(loc='upper right')
    plt.title(f"{label_orbit} in the Neptune–Triton RTBP")
    plt.grid()
    plt.axis('equal')
    plt.show()

def plot_three_body_orbit_3d(sol, mu, label_orbit="RTBP Orbit", label_primary="Neptune", label_secondary="Triton", label_start="Start", color_orbit='green', color_primary='blue', color_secondary='saddlebrown'):
    """
    Plot the 3D RTBP orbit.
    :param sol: (N x 3) solution array with [x, y, z]
    :param mu: Mass ratio (for plotting Neptune and Triton)
    """
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(sol[:, 0], sol[:, 1], sol[:, 2], color=color_orbit, label=label_orbit)
    ax.scatter(-mu, 0, 0, color=color_primary, s=100, label=label_primary)      # Neptune
    ax.scatter(1 - mu, 0, 0, color=color_secondary, s=50, label=label_secondary)  # Triton
    ax.scatter(sol[0, 0], sol[0, 1], sol[0, 2], color=color_orbit, s=30, label=label_start)  # Start

    ax.set_xlabel("x (non-dimensional)")
    ax.set_ylabel("y (non-dimensional)")
    ax.set_zlabel("z (non-dimensional)")
    ax.set_title(f"{label_orbit} in the Neptune–Triton RTBP")
    ax.legend(loc='upper right')
    plt.show()


# Gravitational parameters for Neptune and Triton
mu_r3bp = mass_ratio(M_neptune, M_triton)
#print(f"Mass ratio (mu) for Neptune-Triton RTBP: {mu_r3bp:.6f}")

# Initialize symbolic gradient functions
U_x_func, U_y_func, U_z_func = get_potential_gradient_funcs()

# Initial Conditions (non-dimensional, rotating frame)
x0, y0, z0 = 0.5, 0.0, 0.0
vx0, vy0, vz0 = 0.0, 1.0, 0.0
state0 = [x0, y0, z0, vx0, vy0, vz0]

# Integration parameters
t_span = (0, 100)
t_eval = np.linspace(*t_span, 1000)

# Solve ODE
solution = solve_ivp(r3bp_equations, t_span, state0, t_eval=t_eval, args=(mu_r3bp,), method='RK45')

# Extract solution: shape (N, 3)
sol_xyz = np.vstack((solution.y[0], solution.y[1], solution.y[2])).T  # (N, 3)


# Plot the RTBP orbit
if __name__ == "__main__":
    plot_three_body_orbit_3d(
        sol_xyz,
        mu_r3bp,
        label_orbit="Triton's RTBP Orbit",
        label_primary="Neptune",
        label_secondary="Triton",
        label_start="Initial Position"
    )

