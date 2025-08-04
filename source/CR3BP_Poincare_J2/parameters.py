"""
parameters.py

Input the parameters for the CR3BP_Poincare_J2 to run the simulations.
These parameters will be used to create the equations and values for the range
of Jacobi constants and initial conditions for the x0 sweep.

Author: Blake T. Johnson
Thesis Project
"""
# === Jacobi Constant Sweep ===
C0 = 3.0144       # Initial Jacobi constant
CF = 3.0144      # Final Jacobi constant
dC = 0.1000       # Jacobi constant increment

# === x0 Sweep Range (Non-Dimensional) ===
XI = -0.20       # Initial x0 (non-dimensional)
XF = 1.20       # Final x0 (non-dimensional)
DX = 0.1000       # Step size for x0 (non-dimensional)

# === Minimum safe distances from each body (km) ===
min_distance_neptune_km = 0.0   # Above surface
min_distance_triton_km  = 0.0   # Above surface

# === Initial y0 Value (Non-Dimensional) ===
YI = 0.0         # Default: 0.0 (can set to sqrt(3)/2 for L4/L5

# === Time Parameters (Dimensional) ===
tlim_sec = 20000.0   # Integration time limit in seconds
dt_sec = 10.0      # Integration timestep in seconds


# === Simulation Flags ===
J2_enabled = False  # If True, include J2 perturbation in the equations of motion

# === Mapping Mode ===
# "highres" = event detection, fine DX, exact crossings
# "global"  = fixed step, midpoint interpolation, professor-style
mapping_mode = "highres"


