"""
jacobi_constant_03.py

This script computes the Jacobi Constant in the Circular Restricted Three-Body Problem (CR3BP)
for the Neptune–Triton system. The Jacobi Constant is an energy-like integral of motion and 
helps define forbidden regions for a third body in a rotating frame.

References:
- Murray & Dermott, Solar System Dynamics, Cambridge University Press, 1999, Ch. 3
- Equation: Jacobi Constant (Murray p. 68, Eq. 3.23)
- Effective Potential: (Murray p. 67, Eq. 3.22)

Usage:
    python jacobi_constant_03.py
"""

import numpy as np
import sys

sys.path.append('/Users/blakejohnson/Documents/Thesis/Three Body Problem/neptune_triton_models/r3bp')
from r3bp_calculations_02 import state0, mu_r3bp  # Import initial state and mu

def jacobi_constant(state, mu):
    """
    Calculate the Jacobi constant for a given state vector in the CR3BP.
    
    Parameters:
        state : list or array
            [x, y, z, vx, vy, vz] position and velocity in the rotating frame
        mu : float
            Mass parameter, defined as mu = m2 / (m1 + m2)

    Returns:
        C : float
            The Jacobi constant value
        
    Equation Reference:
        C = 2U - (vx^2 + vy^2 + vz^2)
        where U is the effective potential:
        U = 0.5*(x^2 + y^2 + z^2) + (1 - mu)/r1 + mu/r2
        r1 = distance to primary (Neptune), r2 = distance to secondary (Triton)
        Source: Murray p. 68, Eq. 3.23
    """
    x, y, z, vx, vy, vz = state
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x + mu - 1)**2 + y**2 + z**2)
    U = 0.5 * (x**2 + y**2 + z**2) + (1 - mu) / r1 + mu / r2
    C = 2 * U - (vx**2 + vy**2 + vz**2)
    return C

C = jacobi_constant(state0, mu_r3bp)
print(f"Jacobi constant (C): {C:.6f}")



