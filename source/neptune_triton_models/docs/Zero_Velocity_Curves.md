# Zero-Velocity Curves Documentation

## Theoretical Foundation

Zero-velocity curves (also called zero-velocity surfaces in 3D) are fundamental structures in the Circular Restricted Three-Body Problem that define the boundaries between accessible and forbidden regions of motion. They represent surfaces in configuration space where a test particle must have zero velocity in the rotating reference frame.

## Mathematical Definition

### Basic Concept
At any point where the particle has zero velocity in the rotating frame:
```
v² = ẋ² + ẏ² + ż² = 0
```

From the definition of the Jacobi constant:
```
C = 2U(x,y,z) - v²
```

When v² = 0, we have:
```
C = 2U(x,y,z)
```

Therefore, zero-velocity curves are contours of the effective potential where:
```
2U(x,y,z) = C = constant
```

### Effective Potential
The effective potential in the rotating frame is:
```
U(x,y,z) = ½(x² + y²) + (1-μ)/r₁ + μ/r₂
```

Where:
- `r₁ = √[(x+μ)² + y² + z²]` = distance to Neptune
- `r₂ = √[(x-1+μ)² + y² + z²]` = distance to Triton  
- `μ = M_triton/(M_neptune + M_triton)` ≈ 2.09×10⁻⁴

### Motion Constraints
The zero-velocity condition defines:
- **Accessible regions**: Where C ≤ 2U (motion possible)
- **Forbidden regions**: Where C > 2U (motion impossible)
- **Boundary surfaces**: Where C = 2U (turning points)

## Physical Interpretation

### Energy Perspective
In the rotating frame, total energy per unit mass is:
```
E_rotating = ½v² + U_effective
```

At zero-velocity surfaces:
- **Kinetic energy**: Zero in rotating frame
- **Potential energy**: Maximum for given Jacobi constant
- **Turning points**: Particle direction reverses

### Classical Analogy
Similar to classical mechanics where a particle bounces off potential barriers:
- **Hills in potential**: Create forbidden regions
- **Valleys in potential**: Allow continued motion
- **Saddle points**: Provide gateways between regions

## Critical Jacobi Constant Values

The topology of zero-velocity curves changes dramatically at specific values of the Jacobi constant corresponding to the Lagrange points.

### C₁ (Jacobi Constant at L₁)
- **Value**: C₁ ≈ 3.172 for Neptune-Triton system
- **Topology**: Separate forbidden regions around each primary
- **Motion type**: Particle trapped near one primary only

### C₂ (Jacobi Constant at L₂)  
- **Value**: C₂ ≈ 3.032
- **Topology**: Forbidden regions merge, opening L₁ gateway
- **Motion type**: Transfer between primaries becomes possible

### C₃ (Jacobi Constant at L₃)
- **Value**: C₃ ≈ 3.012  
- **Topology**: Access to exterior regions beyond L₂
- **Motion type**: Escape to larger distances possible

### C₄ = C₅ (Jacobi Constant at L₄, L₅)
- **Value**: C₄ = C₅ ≈ 3.000
- **Topology**: Complete access to triangular point regions
- **Motion type**: Motion around both primaries simultaneously

## Visualization and Computation

### 2D Contour Plots
For motion in the xy-plane (z = 0):

```python
import numpy as np
import matplotlib.pyplot as plt

def plot_zero_velocity_curves(C_values, mu, xlim=(-2, 2), ylim=(-2, 2), resolution=1000):
    """Plot zero-velocity curves for given Jacobi constants"""
    
    x = np.linspace(xlim[0], xlim[1], resolution)
    y = np.linspace(ylim[0], ylim[1], resolution)
    X, Y = np.meshgrid(x, y)
    
    # Calculate effective potential
    r1 = np.sqrt((X + mu)**2 + Y**2)
    r2 = np.sqrt((X - 1 + mu)**2 + Y**2)
    
    # Avoid singularities at primary locations
    r1 = np.maximum(r1, 1e-6)
    r2 = np.maximum(r2, 1e-6)
    
    U = 0.5*(X**2 + Y**2) + (1-mu)/r1 + mu/r2
    
    # Plot contours for each Jacobi constant
    fig, ax = plt.subplots(figsize=(10, 8))
    
    for i, C in enumerate(C_values):
        contours = ax.contour(X, Y, 2*U, levels=[C], 
                             colors=f'C{i}', linewidths=2)
        ax.clabel(contours, inline=True, fontsize=10, fmt=f'C={C:.3f}')
    
    # Mark primary positions
    ax.plot(-mu, 0, 'bo', markersize=8, label='Neptune')
    ax.plot(1-mu, 0, 'ro', markersize=6, label='Triton')
    
    # Mark Lagrange points
    L_points = find_lagrange_points(mu)
    for name, coords in L_points.items():
        ax.plot(coords[0], coords[1], 'k+', markersize=8)
        ax.annotate(name, (coords[0], coords[1]), xytext=(5, 5), 
                   textcoords='offset points')
    
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel('x (normalized)')
    ax.set_ylabel('y (normalized)')
    ax.set_title('Zero-Velocity Curves in Neptune-Triton System')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    
    return fig, ax
```

### 3D Surface Visualization
The effective potential can be visualized as a 3D surface:

```python
def plot_effective_potential_3d(mu, xlim=(-1.5, 1.5), ylim=(-1.5, 1.5)):
    """Plot 3D surface of effective potential"""
    
    x = np.linspace(xlim[0], xlim[1], 100)
    y = np.linspace(ylim[0], ylim[1], 100)
    X, Y = np.meshgrid(x, y)
    
    r1 = np.sqrt((X + mu)**2 + Y**2)
    r2 = np.sqrt((X - 1 + mu)**2 + Y**2)
    
    # Limit potential for visualization (avoid singularities)
    r1 = np.maximum(r1, 0.05)
    r2 = np.maximum(r2, 0.05)
    
    U = 0.5*(X**2 + Y**2) + (1-mu)/r1 + mu/r2
    
    # Cap potential for better visualization
    U = np.minimum(U, 10)
    
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection='3d')
    
    surf = ax.plot_surface(X, Y, U, cmap='viridis', alpha=0.8, 
                          linewidth=0, antialiased=True)
    
    # Add contour lines at base
    ax.contour(X, Y, U, levels=10, offset=np.min(U), colors='gray', alpha=0.5)
    
    ax.set_xlabel('x (normalized)')
    ax.set_ylabel('y (normalized)')  
    ax.set_zlabel('Effective Potential U')
    ax.set_title('Effective Potential Landscape')
    
    fig.colorbar(surf, shrink=0.5, aspect=10)
    
    return fig, ax
```

## Regions of Motion

### Classification by Jacobi Constant

#### High Energy (C < C₄)
- **Complete freedom**: Access to all space regions
- **Unrestricted motion**: Can reach infinity
- **Practical significance**: Escape trajectories, interplanetary transfer

#### Moderate Energy (C₄ < C < C₃)
- **Triangular regions accessible**: Motion around L₄, L₅ possible
- **Still exterior access**: Can reach large distances
- **Applications**: Trojan-type orbits, horseshoe trajectories

#### Intermediate Energy (C₃ < C < C₂)
- **L₃ gateway closed**: No access beyond L₃
- **Interior access maintained**: Can still reach L₁ region
- **Bounded exterior motion**: Limited to vicinity of system

#### Lower Energy (C₂ < C < C₁)  
- **L₂ gateway closed**: No access beyond Triton
- **L₁ channel open**: Transfer between primaries possible
- **Confined motion**: Bound to immediate vicinity of primaries

#### Lowest Energy (C > C₁)
- **Complete separation**: Independent motion around each primary
- **No mass transfer**: L₁ gateway closed
- **Individual Hill spheres**: Distinct gravitational domains

### Topology Evolution
As energy decreases (C increases):
1. **Single connected region** (high energy)
2. **Constriction at L₂** (moderate energy)
3. **Constriction at L₃** (intermediate energy)  
4. **Bottleneck at L₁** (lower energy)
5. **Complete separation** (lowest energy)

## Applications in Orbital Mechanics

### Mission Design

#### Transfer Trajectory Planning
- **Energy requirements**: Determine minimum ΔV for transfers
- **Gateway identification**: Locate low-energy transfer routes
- **Trajectory constraints**: Understand accessible regions

#### Station-Keeping
- **Forbidden regions**: Natural boundaries for spacecraft motion
- **Equilibrium neighborhoods**: Regions near Lagrange points
- **Manifold utilization**: Exploit natural dynamics for fuel savings

### Natural Dynamics

#### Asteroid Motion
- **Classification**: Determine long-term behavior based on energy
- **Capture mechanics**: Understand temporary satellite phenomena
- **Resonance effects**: Map stability boundaries

#### Tidal Interactions
- **Mass transfer**: Material flow between binary components
- **Roche lobe**: Limiting surface for mass retention
- **Stream dynamics**: Evolution of transferred material

### Neptune-Triton Applications

#### Debris Analysis
- **Particle fate**: Determine long-term evolution of small debris
- **Capture regions**: Identify zones of temporary capture
- **Escape mechanisms**: Understand loss processes

#### System Evolution  
- **Historical dynamics**: Model past orbital configurations
- **Tidal effects**: Long-term orbital circularization
- **Capture scenarios**: Investigate Triton's acquisition

## Computational Considerations

### Numerical Methods

#### Grid-Based Approach
```python
def compute_zero_velocity_grid(mu, C, xlim, ylim, resolution=500):
    """Compute zero-velocity regions on a grid"""
    
    x = np.linspace(xlim[0], xlim[1], resolution)
    y = np.linspace(ylim[0], ylim[1], resolution)
    X, Y = np.meshgrid(x, y)
    
    r1 = np.sqrt((X + mu)**2 + Y**2)
    r2 = np.sqrt((X - 1 + mu)**2 + Y**2)
    
    U = 0.5*(X**2 + Y**2) + (1-mu)/r1 + mu/r2
    
    # Accessible regions where 2U <= C
    accessible = (2*U <= C)
    forbidden = (2*U > C)
    
    return X, Y, accessible, forbidden, U
```

#### Contour Following
For high-precision boundaries:
```python
def trace_zero_velocity_contour(mu, C, starting_point, step_size=0.001):
    """Trace zero-velocity contour with high precision"""
    
    def potential_diff(point):
        x, y = point
        r1 = np.sqrt((x + mu)**2 + y**2)
        r2 = np.sqrt((x - 1 + mu)**2 + y**2)
        U = 0.5*(x**2 + y**2) + (1-mu)/r1 + mu/r2
        return 2*U - C
    
    # Use contour following algorithm
    # Implementation would use gradient information to follow C = 2U
    pass
```

### Accuracy Considerations
- **Grid resolution**: Balance detail vs. computational cost
- **Singularity handling**: Special treatment near primaries
- **Contour interpolation**: Smooth boundary representation

### Validation Methods
```python
def validate_zero_velocity_curves(X, Y, U, C, tolerance=1e-10):
    """Validate computed zero-velocity curves"""
    
    # Check that boundaries satisfy C = 2U within tolerance
    boundary_points = find_boundary_points(2*U, C)
    errors = []
    
    for point in boundary_points:
        x, y = point
        r1 = np.sqrt((x + mu)**2 + y**2)
        r2 = np.sqrt((x - 1 + mu)**2 + y**2)
        U_point = 0.5*(x**2 + y**2) + (1-mu)/r1 + mu/r2
        error = abs(2*U_point - C)
        errors.append(error)
    
    max_error = max(errors)
    print(f"Maximum boundary error: {max_error:.2e}")
    
    return max_error < tolerance
```

## Advanced Topics

### Hill's Surfaces
Extension to elliptical restricted problem:
- **Time-dependent boundaries**: Varying with orbital phase
- **Resonance effects**: Modified by eccentricity
- **Stability implications**: Changed escape conditions

### Multi-Body Generalizations
Beyond two-primary systems:
- **Additional primaries**: More complex forbidden regions
- **Hierarchical systems**: Nested Hill spheres
- **Galaxy dynamics**: Stellar motion boundaries

### Relativistic Effects
Post-Newtonian corrections:
- **Modified potential**: Additional relativistic terms
- **Precession effects**: Orbital plane rotation
- **Light-speed constraints**: Modified zero-velocity surfaces

## Historical Development

### Classical Period
- **1878**: G.W. Hill introduces concept for lunar theory
- **1892**: H. Poincaré develops mathematical foundation
- **1906**: K.F. Sundman proves convergence of series solutions

### Modern Applications
- **1960s**: Space mission applications
- **1970s**: Computer visualization capabilities
- **1990s**: Manifold theory integration
- **2000s**: Low-energy trajectory design

## References

### Foundational Works
- Hill, G.W. (1878). "Researches in the lunar theory"
- Poincaré, H. (1892). *Les méthodes nouvelles de la mécanique céleste*
- Szebehely, V. (1967). *Theory of Orbits: The Restricted Problem of Three Bodies*

### Modern Treatments
- Murray, C.D. & Dermott, S.F. (1999). *Solar System Dynamics*
- Gómez, G. et al. (2001). *Dynamics and Mission Design Near Libration Points*
- Koon, W.S. et al. (2000). "Low energy transfer to the Moon"

### Applications
- Ross, S.D. (2006). "Statistical theory of interior-exterior transition"
- Anderson, R.L. & Parker, J.S. (2012). "Survey of ballistic transfers to the lunar surface"

---
*Comprehensive guide to forbidden regions and accessible motion in three-body dynamics*  
*Essential tool for trajectory analysis and mission planning*
