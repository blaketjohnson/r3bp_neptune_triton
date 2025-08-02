"""
parameters.py

Input the parameters for the CR3BP_Poincare_J2 to run the simulations.
These parameters will be used to create the equations and values for the range
of Jacobi constants and initial conditions for the x0 sweep.

Author: Blake T. Johnson
Summer Thesis Project, Phase 1
"""

# === Jacobi Constant Sweep ===
C0 = 3.148       # Initial Jacobi constant
CF = 3.148      # Final Jacobi constant
dC = 0.005       # Jacobi constant increment

# === x0 Sweep Range (Non-Dimensional) ===
XI = -0.2       # Initial x0 (non-dimensional)
XF = 1.2       # Final x0 (non-dimensional)
DX = 0.100       # Step size for x0 (non-dimensional)

# === Initial y0 Value (Non-Dimensional) ===
YI = 0.0         # Default: 0.0 (can set to sqrt(3)/2 for L4/L5

# === Time Parameters (Dimensional) ===
tlim_sec = 20000.0   # Integration time limit in seconds
dt_sec = 10.000      # Integration timestep in seconds

# === Simulation Flags ===
J2_enabled = False  # If True, include J2 perturbation in the equations of motion
plot = True # If True, generate plots of the Poincaré sections
combined_plot = True # If True, generate combined plots of the Poincaré sections
use_SOI_as_escape = False  # If True, use Neptune's SOI as escape radius instead of Hill radius

