"""
Restricted Two-Body Problem (R2BP) for the Neptune-Triton System

This module calculates the motion of Neptune and Triton in the R2BP framework.
It uses numerical methods to solve the equations of motion derived from the gravitational interactions.

Author: Blake T. Johnson
Summer Thesis Project, Phase 1
References:
- Curtis, Orbital Mechanics for Engineering Students, Chapter 3
- Murray & Dermott, Solar System Dynamics, p. 66–67
"""

import numpy as np
from scipy.optimize import fsolve
from scipy.constants import G
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from functools import partial
from constants import *

def gravitational_parameter(M1, M2):
    """
    Calculate the gravitational parameter of a two-body system.
    The gravitational parameter (mu) is defined as the product of the gravitational constant (G) and the total mass of the two bodies.
    Curtis p. 63 (2.21)

    :param M1: Mass of the first body (kg)
    :param M2: Mass of the second body (kg)
    :return: Gravitational parameter (m^3/s^2)
    """
    G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
    return G * (M1 + M2)

def radius_at_perigee(a, e):
    """
    Calculate the radius at periapsis (closest point in orbit) for an elliptical orbit.
    Curtis p. 82 (2.73)

    :param a: Semi-major axis of the orbit (m)
    :param e: Eccentricity of the orbit
    :return: Radius at periapsis (m)
    """
    return a * (1 - e)

def semi_minor_axis(a, e):
    """
    Calculate the semi-minor axis of an elliptical orbit.
    Curtis p. 82 (2.76)

    :param a: Semi-major axis of the orbit (m)
    :param e: Eccentricity of the orbit
    :return: Semi-minor axis (m)
    """
    return a * np.sqrt(1 - e**2)

def angular_momentum(mu, r, e):
    """
    Calculate the specific angular momentum of an orbiting body.
    Curtis p. 73 (2.50)

    :param mu: Gravitational parameter of the two-body system (m^3/s^2)
    :param r: Radius at periapsis or apoapsis (m)
    :param e: Eccentricity of the orbit
    :return: Specific angular momentum (kg m^2/s)
    """
    return np.sqrt(mu * r * (1 - e**2))

def velocity_at_perigee(h, r):
    """
    Calculate the velocity at periapsis using specific angular momentum.
    Curtis p. 69 (2.31)

    :param h: Specific angular momentum (kg m^2/s)
    :param r: Radius at periapsis (m)
    :return: Velocity at periapsis (m/s)
    """
    return h / r

def specific_orbital_energy(mu, a):
    """
    Calculate the specific orbital energy of an orbiting body.
    Curtis p. 83 (2.80)

    :param mu: Gravitational parameter of the two-body system (m^3/s^2)
    :param a: Semi-major axis of the orbit (m)
    :return: Specific orbital energy (J/kg)
    """
    return -mu / (2 * a)

def orbital_period(a, mu):
    """
    Calculate the orbital period using Kepler's Third Law.
    Curtis p. 84 (2.83)

    :param a: Semi-major axis of the orbit (m)
    :param mu: Gravitational parameter of the two-body system (m^3/s^2)
    :return: Orbital period (s)
    """
    return (2 * np.pi / np.sqrt(mu)) * (a ** (3 / 2))

def initial_conditions_perigee(rp, vp, M1, M2):
    '''
    Calculate the initial conditions for a two-body problem at perigee.
    :param rp: Radius at perigee (m)
    :param vp: Velocity at perigee (m/s)
    :param M1: Mass of the first body (kg)
    :param M2: Mass of the second body (kg)
    :return: Initial conditions [x0, y0, z0, vx0, vy0, vz0]
    '''
    r0 = np.array([rp * 1000, 0, 0])  # Position (m)
    v0 = np.array([0, vp, 0])  # Velocity (m/s)
    f0 = [r0[0], r0[1], r0[2], v0[0], v0[1], v0[2]]
    return f0

def two_body_equation_of_motion(f, t, mu):
    '''
    Define the equations of motion for a two-body problem.
    :param initial_conditions: Initial conditions [x0, y0, z0, vx0, vy0, vz0]
    :param t: Time (not used in this function, but required by odeint)
    :param mu: Standard gravitational parameter (m^3/s^2)
    :return: Derivatives of the state vector
    '''
    # Unpack initial conditions
    x, y, z, vx, vy, vz = f
    r = np.sqrt(x**2 + y**2 + z**2)
    return [vx, vy, vz, (-mu / r**3) * x, (-mu / r**3) * y, (-mu / r**3) * z]

# Define the two-body equations of motion
#def dfdt(f, t):
    x, y, z, vx, vy, vz = f
    r = np.sqrt(x**2 + y**2 + z**2)
    return [vx, vy, vz, (-mu_2d / r**3) * x, (-mu_2d / r**3) * y, (-mu_2d / r**3) * z]

def solve_two_body_problem(orbital_period, time_step, t0, two_body_equation_of_motion, initial_conditions):
    '''
    Solve the two-body problem using numerical integration.
    :param orbital_period: Orbital period of the satellite (s)
    :param time_step: Time step for integration (s)
    :param t0: Initial time (s)
    :param two_body_equation_of_motion: Function defining the equations of motion
    :param initial_conditions: Initial conditions [x0, y0, z0, vx0, vy0, vz0]
    :return: (solution array in km, time array)
    '''
    tf = orbital_period * 10
    dt = time_step
    tspan = np.arange(t0, tf, dt)
    from scipy.integrate import odeint
    sol = odeint(two_body_equation_of_motion, initial_conditions, tspan)
    sol_km = sol[:, :3] / 1000  # Convert to km
    return sol_km, tspan

def plot_two_body_orbit(sol_km, label_orbit="Orbit", label_primary="Primary", label_start="Start", color_orbit='saddlebrown', color_primary='blue'):
    """
    Plot the 2D orbit from the integrated solution.
    :param sol_km: Solution array (N x 3) with positions in km
    :param label_orbit: Label for the orbit line
    :param label_primary: Label for the primary body
    :param label_start: Label for the starting point
    :param color_orbit: Color for the orbit line
    :param color_primary: Color for the primary body
    """
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(sol_km[:, 0], sol_km[:, 1], color=color_orbit, label=label_orbit)
    plt.scatter(0, 0, color=color_primary, s=50, label=label_primary)
    plt.scatter(sol_km[0, 0], sol_km[0, 1], color=color_orbit, s=30, label=label_start)
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.legend(loc='upper right')
    plt.title(f"{label_orbit} around {label_primary}")
    plt.grid()
    plt.axis('equal')
    plt.show()

# Gravitational parameter for Neptune-Triton system
mu_r2bp = gravitational_parameter(M_neptune, M_triton) 

# Compute the radius at perigee in kilometers and meters
radius_perigee_r2bp = radius_at_perigee(a_triton_km, e_triton)  # km
radius_perigee_r2bp_meters = radius_perigee_r2bp * 1000  # m

# Compute the semi-minor axis in km
semi_minor_axis_r2bp = semi_minor_axis(a_triton_km, e_triton)  # km

# Compute the specific angular momentum at perigee
angular_momentum_r2bp = angular_momentum(mu_r2bp, radius_perigee_r2bp_meters, e_triton)  # (m^2/s)

# Compute the velocity at perigee
velocity_triton_perigee = velocity_at_perigee(angular_momentum_r2bp, radius_perigee_r2bp_meters)  # (m/s)

# Compute the specific orbital energy
specific_orbital_energy_r2bp = specific_orbital_energy(mu_r2bp, a_triton_meters)  #  (J/kg)

# Orbital period for Triton
T_triton_r2bp = orbital_period(a_triton_meters, mu_r2bp)  # s

# Compute the Initial conditions for Triton at perigee
initial_conditions_r2bp = initial_conditions_perigee(radius_perigee_r2bp, velocity_triton_perigee, M_neptune, M_triton) # (m, m/s)

# Define the two-body equations of motion for the restricted two-body problem (R2BP)
equations_of_motion_r2bp = two_body_equation_of_motion(initial_conditions_r2bp, 100, mu_r2bp)

# Create the ODE function for the two-body problem
ode_func = partial(two_body_equation_of_motion, mu=mu_r2bp)

# Solve the equations of motion for Triton
equations_of_motion_triton, tspan = solve_two_body_problem(T_triton_r2bp, 100, 0, ode_func, initial_conditions_r2bp)


#print the results
print(f"Gravitational parameter (mu) for Neptune-Triton system: {mu_r2bp:.2e} m^3/s^2")
print(f"Radius at periapsis for Triton: {radius_perigee_r2bp:.2f} km")
print(f"Semi-minor axis for Triton: {semi_minor_axis_r2bp:.2f} km")
print(f"Angular momentum at periapsis for Triton: {angular_momentum_r2bp:.2e} kg m^2/s")
print(f"Velocity at periapsis for Triton: {velocity_triton_perigee:.2f} m/s")
print(f"Specific orbital energy for Triton: {specific_orbital_energy_r2bp:.2e} J/kg")
print(f"Orbital period for Triton: {T_triton_r2bp:.2f} s")

#Plot the Restricted Two Boyd orbit of Triton around Neptune
if __name__ == "__main__":
    plot_two_body_orbit(
        equations_of_motion_triton,
        label_orbit="Triton's Orbit",
        label_primary="Neptune",
        label_start="Triton at Perigee"
    )






