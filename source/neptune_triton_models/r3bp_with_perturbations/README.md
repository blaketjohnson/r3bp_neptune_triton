# R3BP with Perturbations

## Overview
This module extends the classical Circular Restricted Three-Body Problem (CR3BP) to include realistic perturbations, primarily Neptune's J₂ oblateness effect. This provides a more accurate model of the Neptune-Triton system dynamics while maintaining the computational advantages of the R3BP framework.

## Key Features
- **J₂ Oblateness Perturbation**: Neptune's quadrupole gravitational field
- **Enhanced Dynamics**: More realistic orbital evolution
- **Stability Analysis**: Perturbed Lagrange point stability
- **Zero-Velocity Curves**: Modified by oblateness effects
- **Comparison Tools**: Direct comparison with classical R3BP results

## Core Modules

### Primary Dynamics
- **`r3bp_j2_01.py`**: Main R3BP+J₂ implementation
  - Equations of motion with oblateness perturbation
  - Enhanced effective potential including J₂ terms
  - Numerical integration with adaptive step sizing
  - Direct comparison with classical R3BP

### Equilibrium Analysis
- **`lagrange_points_j2_02.py`**: Perturbed Lagrange points
  - Modified equilibrium positions due to J₂
  - Iterative correction from classical solutions
  - Validation of perturbation theory accuracy

- **`stability_lagrange_j2_04.py`**: Stability with J₂ perturbation
  - Linearized dynamics around perturbed equilibria
  - Eigenvalue analysis including oblateness effects
  - Comparison of stability characteristics

### Phase Space Structure
- **`zero_velocity_curves_j2_03.py`**: Modified zero-velocity surfaces
  - Jacobi-like constant for perturbed system
  - Forbidden region boundaries with J₂ effects
  - Visualization of energy-dependent accessibility

### Utilities
- **`constants.py`**: Physical constants including J₂ coefficients
  - Neptune oblateness parameters
  - Dimensional and non-dimensional conversions
  - Validation against observational data

- **`non_dimensionalizer.py`**: Unit conversion tools
  - Consistent non-dimensionalization scheme
  - Scaling for J₂ perturbation terms
  - Dimensional analysis utilities

## Physical Parameters

### Neptune Oblateness
- **J₂ coefficient**: 3.411 × 10⁻³
- **Equatorial radius**: 24,764 km  
- **Polar flattening**: ~2%
- **Relative magnitude**: J₂ ≈ 0.3% of central field at Triton's distance

### System Scaling
- **Distance unit**: Neptune-Triton separation (354,759 km)
- **Time unit**: Orbital period/(2π) (0.934 days)
- **Mass unit**: Neptune + Triton total mass
- **Non-dimensional Neptune radius**: R_nd ≈ 0.0698

## Mathematical Framework

### Perturbed Effective Potential
The classical R3BP potential is enhanced with J₂ terms:

```
U = U_classical + U_J2

U_J2 = -J₂(R/r₁)²[1 - 3z²/r₁²]/2r₁
```

Where:
- r₁ = distance from Neptune
- R = Neptune's equatorial radius
- z = coordinate perpendicular to orbital plane

### Equations of Motion
```
ẍ - 2ẏ = ∂U/∂x
ÿ + 2ẋ = ∂U/∂y  
z̈ = ∂U/∂z
```

With gradient terms including both classical and J₂ contributions.

## Usage Examples

### Basic Integration with J₂
```python
from r3bp_j2_01 import equations_of_motion_j2
from constants import J2_neptune, mu_nd
from non_dimensionalizer import R_nd

# Initial conditions
state0 = [x0, y0, z0, vx0, vy0, vz0]

# Integration parameters  
params = (mu_nd, J2_neptune, R_nd)

# Solve with J₂ perturbation
sol = solve_ivp(equations_of_motion_j2, [0, t_final], state0, 
                args=params, rtol=1e-12)
```

### Lagrange Point Comparison
```python
from lagrange_points_j2_02 import find_lagrange_points_j2
from lagrange_points_04 import find_lagrange_points  # classical

# Classical positions
L_classical = find_lagrange_points(mu)

# Perturbed positions
L_j2 = find_lagrange_points_j2(mu, J2_neptune, R_nd)

# Compare shifts
delta_L = L_j2 - L_classical
```

### Zero-Velocity Curve Modification
```python
from zero_velocity_curves_j2_03 import plot_zvc_j2
from zero_velocity_curves_06 import plot_zero_velocity_curves

# Plot comparison
plot_zero_velocity_curves(C, mu)  # classical
plot_zvc_j2(C, mu, J2_neptune, R_nd)  # with J₂
```

## Perturbation Effects

### Magnitude Assessment
- **At Triton's orbit**: J₂ perturbation ≈ 0.3% of central force
- **Lagrange point shifts**: ~10⁻⁴ in normalized units
- **Stability modifications**: Small but measurable eigenvalue changes
- **Long-term evolution**: Cumulative effects can be significant

### Physical Interpretation
J₂ perturbations cause:
- **Precession effects**: Orbital plane and argument changes
- **Resonance shifts**: Modified commensurability conditions  
- **Stability boundaries**: Altered regular/chaotic motion boundaries
- **Transport modification**: Changed phase space flux rates

### Computational Considerations
- **Integration time**: ~20% increase due to additional force terms
- **Accuracy requirements**: Higher precision needed for long-term evolution
- **Convergence**: Iterative methods for equilibrium finding
- **Validation**: Systematic comparison with classical results

## Applications

### Research Questions
1. How does Neptune's oblateness affect Triton's long-term orbital evolution?
2. Do perturbed Lagrange points remain stable for realistic system parameters?
3. What are the observational signatures of J₂ effects in the system?
4. How do perturbations modify chaotic transport in phase space?

### Future Extensions
- **Additional perturbations**: Tidal forces, solar perturbations
- **Higher-order terms**: J₄, J₆ coefficients
- **Triton's oblateness**: Mutual tidal effects
- **Relativistic corrections**: Post-Newtonian dynamics

## References
- Murray, C.D. & Dermott, S.F. (1999). *Solar System Dynamics*, Ch. 4
- Jacobson, R.A. (2009). "The orbits of the Neptunian satellites"
- Lainey, V. et al. (2006). "Natural satellite ephemerides"

---
*Advanced perturbation analysis for Neptune-Triton thesis project*  
*Author: Blake T. Johnson, 2025*
