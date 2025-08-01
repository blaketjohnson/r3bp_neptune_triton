"""
parameters.py

Input the parameters for the CR3BP_Poincare_J2 to run the simulations.
These parameters will be used to create the equations and values for the range
of Jacobi constants and initial conditions for the x0 sweep.

Author: Blake T. Johnson
Summer Thesis Project, Phase 1
"""

# === Jacobi Constant Sweep ===
C0 = 3.02       # Initial Jacobi constant
CF = 3.02       # Final Jacobi constant
dC = 0.001       # Jacobi constant increment

# === x0 Sweep Range (Non-Dimensional) ===
XI = 0.205       # Initial x0 (non-dimensional)
XF = 0.205       # Final x0 (non-dimensional)
DX = 0.005       # Step size for x0 (non-dimensional)

# === Time Parameters (Dimensional) ===
tlim_sec = 10000.0   # Integration time limit in seconds
dt_sec = 0.001      # Integration timestep in seconds

# === Simulation Flags ===
J2_enabled = False  # If True, include J2 perturbation in the equations of motion
plot = True # If True, generate plots of the Poincaré sections
combined_plot = True # If True, generate combined plots of the Poincaré sections
use_SOI_as_escape = False  # If True, use Neptune's SOI as escape radius instead of Hill radius

