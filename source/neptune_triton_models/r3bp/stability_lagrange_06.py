"""
stability_lagrange_06.py

This script analyzes the stability of the Lagrange Points (L1 through L5) in the Circular 
Restricted Three-Body Problem (CR3BP) for the Neptune–Triton system. Stability is 
assessed via eigenvalue analysis of the linearized system matrix near each Lagrange point.

Key Features:
- Computes second partial derivatives of the effective potential.
- Constructs the Jacobian matrix of the linearized equations of motion.
- Calculates eigenvalues to classify each Lagrange point as stable or unstable.

References:
- Murray & Dermott, Solar System Dynamics, Cambridge University Press, 1999.
  - Effective Potential Derivatives: Ch. 3, pp. 67–70
  - Stability Criteria: p. 83, Eq. 3.32
- Uses the Lagrange points computed in lagrange_points_04.py.

Usage:
    python stability_lagrange_06.py
"""


import numpy as np
import math
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
from scipy.optimize import root
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sympy as sp
import sys
from zero_velocity_curves_05 import *


# Function to calculate and return the second partial derivatives Uxx, Uxy, Uyy
def dU():
    x, y, mu = sp.symbols('x y mu')
    r1 = sp.sqrt((x + mu)**2 + y**2)
    r2 = sp.sqrt((x + mu - 1)**2 + y**2)
    U_sym = 0.5 * (x**2 + y**2) + (1 - mu) / r1 + mu / r2
    
    Uxx = sp.diff(U_sym, x, x)
    Uxy = sp.diff(U_sym, x, y)
    Uyy = sp.diff(U_sym, y, y)
    
    return Uxx, Uxy, Uyy

# Function to calculate the second partial derivatives Uxx, Uxy, Uyy with given values
def ddU(x_val, y_val, mu_val):
    Uxx, Uxy, Uyy = dU()  # Calculate the symbolic derivatives
    Uxx_val = Uxx.subs({sp.symbols('x'): x_val, sp.symbols('y'): y_val, sp.symbols('mu'): mu_val})
    Uxy_val = Uxy.subs({sp.symbols('x'): x_val, sp.symbols('y'): y_val, sp.symbols('mu'): mu_val})
    Uyy_val = Uyy.subs({sp.symbols('x'): x_val, sp.symbols('y'): y_val, sp.symbols('mu'): mu_val})
    
    # Convert to real part or keep as complex
    Uxx_val = complex(Uxx_val)
    Uxy_val = complex(Uxy_val)
    Uyy_val = complex(Uyy_val)
    
    return Uxx_val, Uxy_val, Uyy_val

# Function to round small numbers to zero based on a threshold
def round_small_numbers_to_zero(value, threshold=1e-10):
    if abs(value) < threshold:
        return 0
    else:
        return value

# Function to round small components of eigenvalues
def round_eigenvalues(eigenvalues, threshold=1e-10):
    rounded_eigenvalues = []
    for eigenvalue in eigenvalues:
        real_part = round_small_numbers_to_zero(sp.re(eigenvalue), threshold)
        imag_part = round_small_numbers_to_zero(sp.im(eigenvalue), threshold)
        rounded_eigenvalues.append(real_part + imag_part * sp.I)
    return rounded_eigenvalues

# Function to calculate the eigenvalues of the system matrix
def eigenvalues(x_val, y_val, mu_val):
    Uxx_val, Uxy_val, Uyy_val = ddU(x_val, y_val, mu_val)
    
    # Print Uxx, Uxy, Uyy values for inspection
    print(f"Uxx: {Uxx_val}, Uxy: {Uxy_val}, Uyy: {Uyy_val}")
    
    # Define the matrix A
    A = sp.Matrix([
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [Uxx_val, Uxy_val, 0, 2],
        [Uxy_val, Uyy_val, -2, 0]
    ])
    
    # Print matrix A for inspection
    print("Matrix A:")
    sp.pprint(A, use_unicode=True)
    
    # Calculate the eigenvalues of matrix A
    eigenvalues = A.eigenvals()
    
    return eigenvalues

# Function to determine stability based on eigenvalues
def check_stability(eigenvalues):
    # Check the real parts of the eigenvalues
    for eigenvalue in eigenvalues:
        real_part = sp.re(eigenvalue)
        if real_part > 0:
            return "Unstable"
    return "Stable"

# List of Lagrange points (assuming you have these values pre-calculated)
Lagrange_points = [L1, L2, L3, L4, L5]

# Loop over the Lagrange points and compute eigenvalues for each
for i, L in enumerate(Lagrange_points, start=1):
    x_val, y_val = L[0], L[1]
    print(f"Calculating eigenvalues for L{i} at (x = {x_val}, y = {y_val})")
    
    # Call your original eigenvalues function
    eigenvalues_result = eigenvalues(x_val, y_val, mu_r3bp)
    
    # Round small components of eigenvalues
    rounded_eigenvalues = round_eigenvalues(eigenvalues_result)
    
    # Determine stability based on eigenvalues
    stability = check_stability(rounded_eigenvalues)
    
    # Display the rounded eigenvalues
    print(f"Rounded Eigenvalues for L{i}:")
    for eigenvalue in rounded_eigenvalues:
        sp.pprint(eigenvalue, use_unicode=True)
    
    # Display the stability result
    print(f"L{i} is {stability}\n")