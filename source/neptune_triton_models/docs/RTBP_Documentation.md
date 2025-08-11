# Restricted Three-Body Problem (R3BP) Documentation

## Theoretical Framework

The Circular Restricted Three-Body Problem (CR3BP) is a fundamental model in celestial mechanics that describes the motion of a massless test particle under the gravitational influence of two massive bodies (primaries) that orbit their common center of mass in circular orbits.

## Problem Setup

### Physical Configuration
- **Primary bodies**: Neptune (m₁) and Triton (m₂) in circular orbit
- **Test particle**: Massless satellite or debris
- **Assumptions**: 
  - Circular motion of primaries
  - Test particle doesn't affect primary motion
  - No other forces (purely gravitational)

### Coordinate Systems

#### Inertial Frame
- Origin at system barycenter
- Fixed orientation in space
- Newton's laws apply directly

#### Rotating Frame (Synodic)
- Origin at system barycenter  
- x-axis always points from Neptune to Triton
- z-axis perpendicular to orbital plane
- Frame rotates with angular velocity ω

## Mathematical Formulation

### Non-dimensional Units
To simplify analysis, we use normalized coordinates:
- **Distance unit**: Separation between primaries = 1
- **Time unit**: Orbital period/(2π) = 1
- **Mass unit**: Total mass of primaries = 1
- **Velocity unit**: Distance/Time

### Mass Ratio
```
μ = m₂/(m₁ + m₂) = M_triton/(M_neptune + M_triton) ≈ 2.09 × 10⁻⁴
```

### Primary Positions (Rotating Frame)
- **Neptune**: (-μ, 0, 0)
- **Triton**: (1-μ, 0, 0)
- **Barycenter**: (0, 0, 0)

### Effective Potential
The key quantity in the rotating frame is the effective potential:

```
U(x,y,z) = ½(x² + y²) + (1-μ)/r₁ + μ/r₂ + ½μ(1-μ)
```

Where:
- `r₁ = √[(x+μ)² + y² + z²]` = distance to Neptune
- `r₂ = √[(x-1+μ)² + y² + z²]` = distance to Triton

### Equations of Motion
In the rotating frame, the equations of motion are:

```
ẍ - 2ẏ = ∂U/∂x
ÿ + 2ẋ = ∂U/∂y  
z̈ = ∂U/∂z
```

The terms `2ẏ` and `2ẋ` are Coriolis force components due to frame rotation.

### Partial Derivatives
```
∂U/∂x = x - (1-μ)(x+μ)/r₁³ - μ(x-1+μ)/r₂³
∂U/∂y = y - (1-μ)y/r₁³ - μy/r₂³
∂U/∂z = -(1-μ)z/r₁³ - μz/r₂³
```

## Jacobi Constant (Energy Integral)

### Definition
The Jacobi constant C is the only known integral of motion in the CR3BP:

```
C = 2U(x,y,z) - (ẋ² + ẏ² + ż²)
```

This represents twice the effective potential energy minus the kinetic energy.

### Physical Interpretation
- **Energy boundary**: Defines accessible regions in phase space
- **Motion constraint**: C = constant along trajectories
- **Zero-velocity surfaces**: Boundaries where kinetic energy = 0

### Relationship to Classical Energy
In the inertial frame, the classical energy is:
```
E_inertial = ½(v_x² + v_y² + v_z²) - (1-μ)/r₁ - μ/r₂
```

The Jacobi constant relates to this through the coordinate transformation.

## Lagrange Points (Equilibrium Solutions)

### Definition
Lagrange points are equilibrium positions in the rotating frame where ∇U = 0.

### Collinear Points (L₁, L₂, L₃)
Located on the x-axis, found by solving:
```
∂U/∂x = 0, y = z = 0
```

#### L₁ Point (Between Primaries)
- **Location**: Between Neptune and Triton
- **Stability**: Unstable (saddle point)
- **Physical significance**: Gateway for mass transfer

#### L₂ Point (Beyond Triton)  
- **Location**: On Triton side, beyond Triton
- **Stability**: Unstable (saddle point)
- **Applications**: Spacecraft staging area

#### L₃ Point (Beyond Neptune)
- **Location**: On Neptune side, beyond Neptune  
- **Stability**: Unstable (saddle point)
- **Least relevant**: Far from system dynamics

### Triangular Points (L₄, L₅)
- **Location**: Form equilateral triangles with primaries
- **Coordinates**: (½-μ, ±√3/2, 0)
- **Stability**: Linearly stable for μ < μ_critical ≈ 0.0385

For Neptune-Triton system: μ ≈ 2.09×10⁻⁴ << μ_critical, so L₄ and L₅ are stable.

## Numerical Integration Considerations

### Symplectic Integration
The Hamiltonian nature of the CR3BP makes symplectic integrators ideal:
- **Preserve energy** over long integration times
- **Maintain phase space structure**
- **Examples**: Verlet, Ruth-Forest, Yoshida methods

### Adaptive Methods
For high-accuracy applications:
- **Runge-Kutta-Fehlberg (RKF45)**
- **Dormand-Prince (DOPRI)**
- **Variable step-size** for efficiency

### Regularization Techniques
Near primaries, gravitational forces become singular:
- **Levi-Civita regularization**: For binary encounters
- **Kustaanheimo-Stiefel**: Quaternion-based regularization
- **Switching methods**: Different integrators in different regions

## Physical Applications

### Mission Design
- **Halo orbits**: Periodic orbits around L₁, L₂
- **Lissajous trajectories**: Quasi-periodic motion near Lagrange points
- **Weak stability boundary theory**: Low-energy transfers

### Natural Dynamics
- **Trojan asteroids**: Objects at Jupiter's L₄, L₅ points
- **Tidal interactions**: Mass transfer in binary systems
- **Ring dynamics**: Particle motion in planetary rings

### Neptune-Triton Specifics
- **Triton's retrograde orbit**: Captured object dynamics
- **Tidal evolution**: Long-term orbital changes
- **Debris dynamics**: Small particle motion in system

## Stability Analysis

### Linear Stability
Near equilibrium points, linearize equations:
```
δẍ - 2δẏ = U_xx δx + U_xy δy + U_xz δz
δÿ + 2δẋ = U_yx δx + U_yy δy + U_yz δz
δz̈ = U_zx δx + U_zy δy + U_zz δz
```

### Characteristic Equation
The stability is determined by eigenvalues of the linearized system.

#### For Collinear Points
Eigenvalues are typically:
- Two real (one positive, one negative): **Unstable**
- Two pure imaginary: **Center-type behavior**

#### For Triangular Points  
When linearly stable (μ < μ_critical):
- All eigenvalues pure imaginary: **Stable oscillations**
- Frequencies depend on mass ratio

## Poincaré Surface of Section

### Definition
A systematic way to study phase space structure by recording trajectory crossings of a chosen surface.

### Standard Choice
For CR3BP, commonly use:
- **Surface**: y = 0 (x-axis crossing)
- **Direction**: ẏ > 0 (upward crossing)
- **Coordinates**: (x, ẋ) at crossing

### Information Content
- **Fixed points**: Periodic orbits
- **Invariant curves**: Quasi-periodic motion (KAM tori)
- **Chaotic regions**: Stochastic trajectories
- **Separatrices**: Boundaries between motion types

## Connection to Real Systems

### Perturbations
Real systems deviate from CR3BP through:
- **Eccentricity**: Non-circular primary orbit
- **Inclination**: Non-planar motion
- **Oblateness**: J₂, J₄ gravitational harmonics  
- **Solar radiation pressure**: For small particles
- **Tidal forces**: From other bodies

### Neptune-Triton Enhancements
Specific corrections for realism:
- **Neptune's J₂**: Significant oblateness effect
- **Solar perturbations**: Third-body effects
- **Triton's orbital evolution**: Tidal circularization
- **Relativistic corrections**: Post-Newtonian terms

## Computational Implementation

### State Vector
```
state = [x, y, z, ẋ, ẏ, ż]
```

### Equations as First-Order System
```python
def cr3bp_equations(t, state, mu):
    x, y, z, vx, vy, vz = state
    
    r1 = sqrt((x + mu)**2 + y**2 + z**2)
    r2 = sqrt((x - 1 + mu)**2 + y**2 + z**2)
    
    ax = x + 2*vy - (1-mu)*(x+mu)/r1**3 - mu*(x-1+mu)/r2**3
    ay = y - 2*vx - (1-mu)*y/r1**3 - mu*y/r2**3
    az = -(1-mu)*z/r1**3 - mu*z/r2**3
    
    return [vx, vy, vz, ax, ay, az]
```

### Integration Parameters
- **Relative tolerance**: 1e-12 for high accuracy
- **Absolute tolerance**: 1e-15
- **Max step size**: 0.01 (in normalized units)

## References

### Classical Texts
- Szebehely, V. (1967). *Theory of Orbits: The Restricted Problem of Three Bodies*
- Murray, C.D. & Dermott, S.F. (1999). *Solar System Dynamics*
- Moulton, F.R. (1914). *An Introduction to Celestial Mechanics*

### Modern Applications  
- Koon, W.S. et al. (2000). "Low Energy Transfer to the Moon"
- Gómez, G. et al. (2001). *Dynamics and Mission Design Near Libration Points*
- Parker, J.S. & Anderson, R.L. (2013). *Low-Energy Lunar Trajectory Design*

### Neptune-Triton Specific
- McKinnon, W.B. et al. (1995). "Origin and Evolution of Triton"
- Agnor, C.B. & Hamilton, D.P. (2006). "Neptune's capture of its moon Triton"

---
*Comprehensive theoretical foundation for three-body orbital mechanics*  
*Graduate thesis research documentation*
