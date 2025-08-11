# Jacobi Constant Documentation

## Theoretical Foundation

The Jacobi Constant is the only known integral of motion in the Circular Restricted Three-Body Problem (CR3BP), representing a generalized energy that remains constant along any trajectory in the rotating reference frame.

## Mathematical Definition

### Basic Formula
In the rotating (synodic) coordinate system, the Jacobi constant is defined as:

```
C = 2U(x,y,z) - v²
```

Where:
- `U(x,y,z)` = Effective potential in the rotating frame
- `v² = ẋ² + ẏ² + ż²` = Velocity magnitude squared

### Effective Potential
The effective potential combines gravitational and centrifugal effects:

```
U(x,y,z) = ½(x² + y²) + (1-μ)/r₁ + μ/r₂ + ½μ(1-μ)
```

Components:
- `½(x² + y²)`: Centrifugal potential due to frame rotation
- `(1-μ)/r₁`: Gravitational potential from primary 1 (Neptune)
- `μ/r₂`: Gravitational potential from primary 2 (Triton)
- `½μ(1-μ)`: Constant term (often omitted)

### Distance Functions
```
r₁ = √[(x+μ)² + y² + z²]    # Distance to Neptune
r₂ = √[(x-1+μ)² + y² + z²]  # Distance to Triton
```

Where μ = mass ratio = M₂/(M₁+M₂) ≈ 2.09×10⁻⁴ for Neptune-Triton.

## Physical Interpretation

### Energy Conservation
The Jacobi constant represents a modified energy that accounts for:
1. **Kinetic energy** in the rotating frame
2. **Gravitational potential energy** from both primaries  
3. **Centrifugal potential energy** due to frame rotation

### Relationship to Classical Energy
In the inertial frame, the total energy per unit mass is:
```
E = ½v_inertial² - (1-μ)/r₁ - μ/r₂
```

The transformation between frames gives:
```
C = -2E - (1-μ)μ
```

## Zero-Velocity Surfaces

### Definition
Zero-velocity surfaces are regions where the kinetic energy becomes zero:
```
v² = 0  ⟹  C = 2U(x,y,z)
```

At these boundaries, particles must have zero velocity in the rotating frame.

### Physical Significance
- **Forbidden regions**: Areas where C > 2U are inaccessible
- **Accessible regions**: Areas where C ≤ 2U allow motion
- **Tunneling boundaries**: Define possible trajectory types

### Critical Values
Different ranges of C correspond to different motion types:

#### C > C₁ (C₁ ≈ 3.172 for Neptune-Triton)
- **Separate motion**: Particle trapped near one primary
- **Two forbidden regions**: Around each primary separately

#### C₂ < C < C₁ (C₂ ≈ 3.032)  
- **Connected motion**: Particle can move between primaries
- **Single forbidden region**: Interior to both primary regions

#### C₃ < C < C₂ (C₃ ≈ 3.012)
- **Extended motion**: Access to exterior regions
- **L₂ point accessibility**: Beyond Triton

#### C < C₃
- **Unrestricted motion**: Access to all space regions
- **Escape possible**: To infinity

## Computation Methods

### Direct Calculation
For a given state vector [x, y, z, ẋ, ẏ, ż]:

```python
def jacobi_constant(state, mu):
    x, y, z, vx, vy, vz = state
    
    # Distance calculations
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x - 1 + mu)**2 + y**2 + z**2)
    
    # Effective potential
    U = 0.5*(x**2 + y**2) + (1-mu)/r1 + mu/r2
    
    # Velocity squared
    v_squared = vx**2 + vy**2 + vz**2
    
    # Jacobi constant
    C = 2*U - v_squared
    
    return C
```

### Verification During Integration
The Jacobi constant should remain constant during numerical integration:

```python
def verify_jacobi_conservation(solution, mu, tolerance=1e-10):
    C_values = [jacobi_constant(state, mu) for state in solution.y.T]
    C_variation = max(C_values) - min(C_values)
    
    if C_variation < tolerance:
        print(f"Jacobi constant conserved: ΔC = {C_variation:.2e}")
    else:
        print(f"Warning: Jacobi constant drift: ΔC = {C_variation:.2e}")
    
    return C_values
```

## Applications in Orbital Mechanics

### Mission Design
1. **Energy requirements**: Determine minimum ΔV for transfers
2. **Accessible regions**: Identify reachable zones for spacecraft
3. **Trajectory classification**: Categorize orbit types by C value
4. **Stability boundaries**: Define escape and capture conditions

### Natural Dynamics
1. **Asteroid motion**: Classify trajectories in planetary systems
2. **Tidal streams**: Understand debris evolution
3. **Ring dynamics**: Particle motion boundaries
4. **Binary systems**: Mass transfer and stability analysis

### Neptune-Triton System
1. **Debris analysis**: Small particle fate determination
2. **Capture mechanics**: Historical Triton capture scenarios
3. **Stability regions**: Safe zones for hypothetical missions
4. **Resonance effects**: Long-term orbital evolution

## Numerical Considerations

### Accuracy Requirements
- **Integration tolerance**: Use relative tolerance ≤ 1e-12
- **Conservation check**: Monitor C drift during integration
- **Singularity handling**: Special treatment near primaries

### Common Issues
1. **Round-off errors**: Near primary bodies (small r₁ or r₂)
2. **Integration drift**: Accumulation over long times
3. **Initial conditions**: Ensure physical consistency

### Best Practices
```python
# High-precision calculation
def jacobi_constant_precise(state, mu):
    x, y, z, vx, vy, vz = state
    
    # Use higher precision for critical calculations
    r1_sq = (x + mu)**2 + y**2 + z**2
    r2_sq = (x - 1 + mu)**2 + y**2 + z**2
    
    r1 = np.sqrt(r1_sq)
    r2 = np.sqrt(r2_sq)
    
    # Avoid numerical issues near primaries
    if r1 < 1e-8 or r2 < 1e-8:
        print("Warning: Very close to primary body")
    
    U = 0.5*(x**2 + y**2) + (1-mu)/r1 + mu/r2
    v_squared = vx**2 + vy**2 + vz**2
    
    return 2*U - v_squared
```

## Relationship to Other Integrals

### Classical Integrals
In the two-body problem, we have:
1. **Energy**: E = constant
2. **Angular momentum**: L⃗ = constant  
3. **Laplace-Runge-Lenz vector**: A⃗ = constant

The CR3BP breaks most of these symmetries, leaving only the Jacobi constant.

### Restricted Problem Hierarchy
- **R2BP**: 6 integrals (complete integrability)
- **CR3BP**: 1 integral (Jacobi constant only)
- **General 3BP**: 0 additional integrals (non-integrable)

## Lagrange Point Values

The Jacobi constant at each Lagrange point defines critical energy levels:

### Collinear Points
- **C₁**: Jacobi constant at L₁ (≈ 3.172 for Neptune-Triton)
- **C₂**: Jacobi constant at L₂ (≈ 3.032)  
- **C₃**: Jacobi constant at L₃ (≈ 3.012)

### Triangular Points  
- **C₄ = C₅**: Jacobi constant at L₄, L₅ (≈ 3.000)

These values serve as energy thresholds for different types of motion.

## Visualization Techniques

### Contour Plots
Plot zero-velocity curves for different C values:

```python
def plot_zero_velocity_curves(C_values, mu, grid_size=500):
    x = np.linspace(-2, 2, grid_size)
    y = np.linspace(-2, 2, grid_size)
    X, Y = np.meshgrid(x, y)
    
    for C in C_values:
        # Calculate 2U on the grid
        r1 = np.sqrt((X + mu)**2 + Y**2)
        r2 = np.sqrt((X - 1 + mu)**2 + Y**2)
        U = 0.5*(X**2 + Y**2) + (1-mu)/r1 + mu/r2
        
        # Plot C = 2U contour
        plt.contour(X, Y, 2*U, levels=[C], colors='blue')
```

### 3D Surface Plots
Visualize the effective potential landscape:

```python
def plot_effective_potential(mu):
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    x = np.linspace(-2, 2, 100)
    y = np.linspace(-2, 2, 100)
    X, Y = np.meshgrid(x, y)
    
    r1 = np.sqrt((X + mu)**2 + Y**2)
    r2 = np.sqrt((X - 1 + mu)**2 + Y**2)
    
    # Avoid singularities for visualization
    r1 = np.maximum(r1, 0.1)
    r2 = np.maximum(r2, 0.1)
    
    U = 0.5*(X**2 + Y**2) + (1-mu)/r1 + mu/r2
    
    ax.plot_surface(X, Y, U, cmap='viridis', alpha=0.8)
```

## Historical Development

### Key Contributors
- **C.G.J. Jacobi (1836)**: Discovered the integral for restricted problem
- **G.W. Hill (1878)**: Applied to lunar theory
- **E.W. Brown (1896)**: Lunar motion calculations
- **V. Szebehely (1967)**: Modern formulation and applications

### Modern Applications
- **Space mission design**: Apollo program, Genesis mission
- **Asteroid dynamics**: Trojan asteroids, near-Earth objects  
- **Exoplanet systems**: Multi-planet configurations
- **Galactic dynamics**: Stellar motion in galaxy potentials

## References

### Primary Sources
- Jacobi, C.G.J. (1836). "Sur le mouvement d'un point et sur un cas particulier du problème des trois corps"
- Szebehely, V. (1967). *Theory of Orbits: The Restricted Problem of Three Bodies*
- Murray, C.D. & Dermott, S.F. (1999). *Solar System Dynamics*, Chapter 3

### Applications
- Koon, W.S. et al. (2000). "Low Energy Transfer to the Moon"
- Gómez, G. et al. (2001). *Dynamics and Mission Design Near Libration Points*
- Ross, S.D. (2006). "Statistical theory of interior-exterior transition"

---
*Mathematical foundation for energy analysis in three-body dynamics*  
*Essential tool for Neptune-Triton system orbital mechanics*
