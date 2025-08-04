""" 
Constants for Neptune and Triton Models

This module contains constants and known values for Neptune and its moon Triton, as used in various orbital mechanics calculations.
Based on: Howard Curtis, Orbital Mechanics for Engineering Students, Appendix A

Author: Blake T. Johnson
 Thesis Project
"""

# Constants
k = 0.01720209895  # Gaussian Gravitational Constant
c = 2.99792458e8  # Speed of light in m/s
G = 6.672e-11  # Gravitational Constant in m^3 kg^-1 s^-2
AU = 1.495978707e11  # Astronomical unit in meters
day = 23.9345  # Hours (Appendix of Orbital Mechanics)
M_sun = 1.989e30  # Mass of the Sun in kg

# Neptune
M_neptune = 1.0243e26  # Mass in kg
R_neptune = 2.4764e4  # Radius in km
R_neptune_meters = R_neptune * 10**3  # Radius in meters
a_neptune_km = 4.495e9  # Semi-major axis in km
a_neptune_au = 30.06896348  # Semi-major axis in AU
e_neptune = 0.00858587  # Eccentricity
Inclination_neptune = 1.76917  # Inclination in degrees
SOI_neptune = 86600000  # Sphere of influence in km
J2_neptune = 3411e-6 # Neptune's J2 coefficient (unitless), from Murray & Dermott

# Triton
M_triton = 2.15e22  # Mass in kg
R_triton = 1353  # Radius in km
R_trition_meters = R_triton * 10**3  # Radius in meters
a_triton_km = 354760  # Semi-major axis in km
a_triton_meters = a_triton_km * 1000  # Semi-major axis in meters
a_triton_au = a_triton_meters / AU  # Semi-major axis in AU
e_triton = 0.0004  # Eccentricity
Inclination_triton = 156.834  # Inclination in degrees
T_triton_days = 5.876854  # Orbital period in days (positive value)
T_triton_seconds = T_triton_days * 24 * 3600  # Orbital period in seconds

