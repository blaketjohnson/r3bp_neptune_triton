# Neptune-Triton R3BP Core Models

## Overview
This directory contains the core implementation of Restricted Three-Body Problem (R3BP) dynamics for the Neptune-Triton system. These modules provide fundamental calculations for orbital mechanics, Lagrange points, stability analysis, and phase space visualization.

## Core Modules

### Fundamental Calculations
- **`r2bp_calculations_01.py`**: Two-body problem implementation
  - Basic orbital mechanics for Neptune-Triton system
  - Foundation for R3BP analysis
  - Validation against analytical solutions

- **`r3bp_calculations_02.py`**: Three-body problem dynamics
  - Equations of motion in rotating frame
  - Symbolic gradient calculations for effective potential
  - Numerical integration using scipy.solve_ivp

### Physical Constants
- **`constants.py`**: System parameters and unit conversions
  - Neptune and Triton masses, orbital elements
  - Normalized units for R3BP formulation
  - Physical constants in SI and normalized systems

### Equilibrium and Stability
- **`lagrange_points_04.py`**: L1-L5 Lagrange point calculation
  - Numerical root finding for equilibrium positions
  - Validation against theoretical predictions

- **`lagrange_points_05.py`**: Extended Lagrange point analysis
  - Enhanced precision calculations
  - Comparison of different numerical methods

- **`stability_lagrange_07.py`**: Linear stability analysis
  - Eigenvalue analysis at Lagrange points
  - Classification of stability types (stable/unstable/center)
  - Lyapunov orbits and characteristic frequencies

### Energy and Phase Space
- **`jacobi_contant_03.py`**: Jacobi constant calculations
  - Energy integral in rotating frame
  - Forbidden region boundaries
  - Conservation verification during integration

- **`zero_velocity_curves_06.py`**: Zero-velocity surface visualization
  - Contour plots of effective potential
  - Accessible/forbidden region boundaries
  - Energy-dependent phase space structure

### Poincaré Analysis
- **`poincare_08.py`**: Poincaré surface-of-section analysis
  - Phase space mapping and visualization
  - Periodic orbit detection and classification
  - Chaos indicators and regular motion boundaries

## Usage Examples

### Basic R3BP Integration
```python
from r3bp_calculations_02 import equations_of_motion, mass_ratio
from constants import *

# Set up initial conditions
mu = mass_ratio(M_neptune, M_triton)
initial_state = [x0, y0, z0, vx0, vy0, vz0]

# Integrate trajectory
sol = solve_ivp(equations_of_motion, [0, t_final], initial_state, 
                args=(mu,), dense_output=True, rtol=1e-12)
```

### Lagrange Point Calculation
```python
from lagrange_points_04 import find_lagrange_points
from constants import mu_neptune_triton

L_points = find_lagrange_points(mu_neptune_triton)
print(f"L1 position: {L_points[0]}")
```

### Stability Analysis
```python
from stability_lagrange_07 import linear_stability_analysis

eigenvalues, eigenvectors = linear_stability_analysis(L_point_position, mu)
stability_type = classify_stability(eigenvalues)
```

### Zero-Velocity Curves
```python
from zero_velocity_curves_06 import plot_zero_velocity_curves
from jacobi_contant_03 import jacobi_constant

C = jacobi_constant(state, mu)
plot_zero_velocity_curves(C, mu)
```

## Coordinate Systems

### Normalized Units
- Distance: Neptune-Triton separation = 1
- Time: Orbital period/(2π) = 1  
- Mass: Neptune + Triton = 1
- Velocity: Distance/Time

### Rotating Frame
- Origin at system barycenter
- x-axis through Neptune-Triton line
- z-axis perpendicular to orbital plane
- Frame rotates with orbital frequency ω

### Body Positions
- Neptune: (-μ, 0, 0)
- Triton: (1-μ, 0, 0)
- μ = M_triton/(M_neptune + M_triton) ≈ 2.09×10⁻⁴

## Physical Parameters
- **Mass ratio**: μ ≈ 2.09×10⁻⁴
- **Orbital period**: ~5.88 days
- **Semi-major axis**: ~354,759 km
- **Eccentricity**: ~0.000016 (nearly circular)

## Mathematical Framework

The R3BP equations of motion in the rotating frame:
```
ẍ - 2ẏ = Ωₓ
ÿ + 2ẋ = Ωᵧ  
z̈ = Ωᵣ
```

Where Ω is the effective potential:
```
Ω = ½(x² + y²) + (1-μ)/r₁ + μ/r₂
```

The Jacobi constant (energy integral):
```
C = 2Ω - (ẋ² + ẏ² + ż²)
```

## References
- Curtis, H. (2014). *Orbital Mechanics for Engineering Students*
- Murray, C.D. & Dermott, S.F. (1999). *Solar System Dynamics*
- Szebehely, V. (1967). *Theory of Orbits*

---
*Neptune-Triton thesis research project*  
*Author: Blake T. Johnson, 2025*
