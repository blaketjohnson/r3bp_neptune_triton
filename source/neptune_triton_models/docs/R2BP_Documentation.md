# Restricted Two-Body Problem (R2BP) Documentation

## Theoretical Background

The Restricted Two-Body Problem forms the foundation for understanding orbital mechanics and serves as the starting point for more complex multi-body dynamics. In this formulation, we consider the motion of a massless test particle (satellite) under the gravitational influence of two massive bodies.

## Mathematical Formulation

### Gravitational Force Law
The fundamental governing equation is Newton's law of universal gravitation:

```
F = G * m₁ * m₂ / r²
```

Where:
- `G` = Universal gravitational constant (6.674 × 10⁻¹¹ m³/kg/s²)
- `m₁, m₂` = Masses of the two bodies
- `r` = Distance between the bodies

### Equations of Motion
For a test particle of negligible mass in the gravitational field of two massive bodies, the equations of motion in an inertial frame are:

```
r̈ = -μ₁(r - r₁)/|r - r₁|³ - μ₂(r - r₂)/|r - r₂|³
```

Where:
- `μᵢ = G·mᵢ` = Gravitational parameter
- `rᵢ` = Position vector of body i
- `r` = Position vector of test particle

### Key Orbital Parameters

#### Semi-major Axis (a)
Defines the size of the elliptical orbit:
```
a = -μ/(2E)
```
Where `E` is the specific orbital energy.

#### Eccentricity (e)
Describes the shape of the orbit:
- `e = 0`: Circular orbit
- `0 < e < 1`: Elliptical orbit  
- `e = 1`: Parabolic trajectory
- `e > 1`: Hyperbolic trajectory

#### Specific Angular Momentum (h)
```
h = r × v = √(μ·a·(1-e²))
```

#### Orbital Period (T)
For elliptical orbits (Kepler's Third Law):
```
T = 2π√(a³/μ)
```

#### Specific Orbital Energy (E)
```
E = v²/2 - μ/r = -μ/(2a)
```

## Neptune-Triton System Parameters

### Physical Constants
- **Neptune mass**: 1.024 × 10²⁶ kg
- **Triton mass**: 2.147 × 10²² kg  
- **Semi-major axis**: 354,759 km
- **Orbital period**: 5.877 days
- **Eccentricity**: 0.000016 (nearly circular)

### Derived Parameters
```python
# Gravitational parameter
μ = G * (M_neptune + M_triton) = 6.836 × 10¹⁵ m³/s²

# Angular frequency
n = √(μ/a³) = 1.237 × 10⁻⁵ rad/s

# Orbital velocity
v = √(μ/a) = 4.386 km/s
```

## Computational Implementation

### Numerical Integration
The R2BP equations are typically solved using numerical integration methods such as:

1. **Runge-Kutta Methods**: 4th-order RK provides good balance of accuracy and efficiency
2. **Adaptive Step-Size Methods**: RKF45, Dormand-Prince for automatic error control
3. **Symplectic Integrators**: Preserve energy and angular momentum over long integrations

### Initial Conditions
For the Neptune-Triton system, typical initial conditions assume Triton starts at:
- **Position**: Perigee (closest approach)
- **Velocity**: Perpendicular to position vector
- **Coordinates**: 2D motion in orbital plane

### Validation Methods
Computational results should be validated against:
1. **Energy conservation**: Total energy should remain constant
2. **Angular momentum conservation**: For central force problems
3. **Analytical solutions**: Comparison with Keplerian orbits
4. **Observational data**: Match known orbital parameters

## Applications and Limitations

### Applications
- Satellite trajectory planning
- Planetary mission design  
- Asteroid and comet orbit determination
- Foundation for multi-body problems

### Limitations
- Ignores gravitational perturbations from other bodies
- Assumes point masses (no extended body effects)
- Neglects relativistic effects
- No atmospheric drag or other non-gravitational forces

## Connection to R3BP
The R2BP naturally extends to the Restricted Three-Body Problem by:
1. Adding a third massive body
2. Transforming to a rotating reference frame
3. Including centrifugal and Coriolis forces
4. Modifying the effective potential

The R2BP solution provides:
- Initial orbit estimates for R3BP analysis
- Validation benchmarks for numerical methods
- Physical intuition for multi-body dynamics

## References
- Curtis, H.D. (2014). *Orbital Mechanics for Engineering Students*, 3rd Edition
- Vallado, D.A. (2013). *Fundamentals of Astrodynamics and Applications*
- Prussing, J.E. & Conway, B.A. (2012). *Orbital Mechanics*
- Battin, R.H. (1999). *An Introduction to the Mathematics and Methods of Astrodynamics*

## Mathematical Appendix

### Vector Formulation
In Cartesian coordinates:
```
ẍ = -μx/r³
ÿ = -μy/r³  
z̈ = -μz/r³
```

Where `r = √(x² + y² + z²)`

### Energy Integral
The total energy per unit mass:
```
E = (ẋ² + ẏ² + ż²)/2 - μ/r = constant
```

### Angular Momentum Vector
```
h⃗ = r⃗ × v⃗ = constant vector
```

This conservation leads to motion confined to a plane perpendicular to h⃗.

---
*Theoretical foundation for Neptune-Triton orbital mechanics analysis*  
*Part of graduate thesis research project*
