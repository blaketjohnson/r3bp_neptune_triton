"""
poincare_07.py

This script prepares input values for constructing a Poincaré surface of section 
in the Circular Restricted Three-Body Problem (CR3BP) for the Neptune–Triton system.
Key quantities such as the mass ratio (μ), system dimensions, and the Jacobi Constant 
are printed for use in downstream Poincaré map computations.

Key Features:
- Calculates mass parameter μ for the Neptune–Triton system.
- Computes system characteristic distances including SOI and body radii.
- Outputs the Jacobi constant based on previously defined state vector.

References:
- Murray & Dermott, Solar System Dynamics, Cambridge University Press, 1999.
  - Mass Ratio Definition: Ch. 3, Eq. 3.8
  - Jacobi Constant: p. 68, Eq. 3.23
- Initial conditions and potential function defined in:
  - r2bp_calculations_01.py
  - jacobi_contant_03.py
  - constants.py

Usage:
    python poincare_06.py
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.optimize import root
from constants import *
import sys

# Import necessary parameters and functions
sys.path.append('/Users/blakejohnson/Documents/rtbp_neptune_triton/neptune_triton_models')
from r3bp_calculations_02 import mu_r3bp, gradient_U  # updated for modular naming
from r2bp_calculations_01 import radius_perigee_r2bp
from jacobi_contant_03 import C
distance_neptune_triton = R_neptune + R_triton + radius_perigee_r2bp

print(f"r_max (SOI): {SOI_neptune} km")
print(f"r1_min (Neptune): {R_neptune} km")
print(f"r2_min (Triton): {R_triton} km")
print(f"mu_Neptune_triton: {mu_r3bp}")
print(f"The Jacobi constant is: {C}")
print(f"Distance between Neptune and Triton: {distance_neptune_triton:.2f} km")



