# Stability Analysis Documentation

## Theoretical Framework

Stability analysis in the Circular Restricted Three-Body Problem determines the long-term behavior of trajectories near equilibrium points (Lagrange points) and periodic orbits. This analysis is crucial for understanding which regions of phase space support regular motion versus chaotic dynamics.

## Linear Stability Theory

### Mathematical Foundation
Near an equilibrium point, we linearize the equations of motion by considering small perturbations:

```
δr⃗ = r⃗ - r⃗₀
```

Where r⃗₀ is the equilibrium position.

### Linearized Equations of Motion
The linearized system in the rotating frame becomes:

```
δẍ - 2δẏ = U_xx δx + U_xy δy + U_xz δz
δÿ + 2δẋ = U_yx δx + U_yy δy + U_yz δz
δz̈ = U_zx δx + U_zy δy + U_zz δz
```

Where U_ij represents the second partial derivatives of the effective potential evaluated at the equilibrium point.

### State Space Representation
Converting to first-order form with state vector ζ = [δx, δy, δz, δẋ, δẏ, δż]:

```
ζ̇ = Aζ
```

Where A is the 6×6 characteristic matrix:

```
A = [0   0   0   1   0   0 ]
    [0   0   0   0   1   0 ]
    [0   0   0   0   0   1 ]
    [U_xx U_xy U_xz 0   2   0 ]
    [U_yx U_yy U_yz -2  0   0 ]
    [U_zx U_zy U_zz 0   0   0 ]
```

## Second Derivatives of Effective Potential

### General Expressions
The effective potential is:
```
U(x,y,z) = ½(x² + y²) + (1-μ)/r₁ + μ/r₂
```

The second derivatives are:
```
U_xx = 1 - (1-μ)/r₁³ + 3(1-μ)(x+μ)²/r₁⁵ - μ/r₂³ + 3μ(x-1+μ)²/r₂⁵
U_yy = 1 - (1-μ)/r₁³ + 3(1-μ)y²/r₁⁵ - μ/r₂³ + 3μy²/r₂⁵
U_zz = -(1-μ)/r₁³ + 3(1-μ)z²/r₁⁵ - μ/r₂³ + 3μz²/r₂⁵
U_xy = 3(1-μ)(x+μ)y/r₁⁵ + 3μ(x-1+μ)y/r₂⁵
U_xz = 3(1-μ)(x+μ)z/r₁⁵ + 3μ(x-1+μ)z/r₂⁵
U_yz = 3(1-μ)yz/r₁⁵ + 3μyz/r₂⁵
```

### At Lagrange Points
For each Lagrange point, these derivatives take specific values that determine the local stability properties.

## Stability at Collinear Points (L₁, L₂, L₃)

### Characteristic Equation
At collinear points (y = z = 0), the problem decouples into:
- **In-plane motion**: (x,y) coupled by Coriolis forces
- **Out-of-plane motion**: z uncoupled

### In-Plane Dynamics
The characteristic polynomial for the 4×4 in-plane system is:
```
λ⁴ + λ² + U_xx*U_yy = 0
```

Since U_xx > 0 and U_yy < 0 at collinear points, this gives:
- **Two real eigenvalues**: ±σ (hyperbolic, unstable directions)
- **Two pure imaginary**: ±iω (elliptic, oscillatory directions)

### Out-of-Plane Motion
The z-motion gives:
```
λ² + U_zz = 0
```

Since U_zz < 0 at collinear points:
- **Two pure imaginary eigenvalues**: ±i√(-U_zz)

### Stability Classification
**All collinear points are unstable** with signature:
- **Saddle × Center × Center**: One unstable direction, two stable oscillatory modes

### Numerical Implementation
```python
def stability_collinear_point(L_point, mu):
    """Analyze stability at a collinear Lagrange point"""
    x, y, z = L_point  # y = z = 0 for collinear points
    
    # Distance calculations
    r1 = abs(x + mu)
    r2 = abs(x - 1 + mu)
    
    # Second derivatives at the point
    U_xx = 1 - (1-mu)/r1**3 - mu/r2**3 + 3*(1-mu)*(x+mu)**2/r1**5 + 3*mu*(x-1+mu)**2/r2**5
    U_yy = 1 - (1-mu)/r1**3 - mu/r2**3
    U_zz = -(1-mu)/r1**3 - mu/r2**3
    
    # Characteristic matrix (4x4 for in-plane motion)
    A_inplane = np.array([
        [0, 0, 1, 0],
        [0, 0, 0, 1],  
        [U_xx, 0, 0, 2],
        [0, U_yy, -2, 0]
    ])
    
    # Eigenvalue analysis
    eigenvals = np.linalg.eigvals(A_inplane)
    
    # Out-of-plane frequency
    omega_z = np.sqrt(-U_zz)
    
    return eigenvals, omega_z, U_xx, U_yy, U_zz
```

## Stability at Triangular Points (L₄, L₅)

### Critical Mass Ratio
The triangular points are **linearly stable** if and only if:
```
μ < μ_critical = ½(1 - √(23/27)) ≈ 0.0385
```

For Neptune-Triton: μ ≈ 2.09×10⁻⁴ << μ_critical, so L₄ and L₅ are stable.

### Frequency Analysis
When stable, the motion consists of small-amplitude oscillations with two fundamental frequencies:

#### Short-Period Libration
```
ω₁² = ½[4 - 27μ + √((4-27μ)² - 4(1-μ)²)]
```

For small μ: ω₁ ≈ √(3/4) ≈ 0.866

#### Long-Period Libration  
```
ω₂² = ½[4 - 27μ - √((4-27μ)² - 4(1-μ)²)]
```

For small μ: ω₂ ≈ √(9μ/4) ≈ 1.5√μ

### Physical Interpretation
- **Short period**: ≈ 7.3 time units (≈ 43 days for Neptune-Triton)
- **Long period**: ≈ 10.5/√μ time units (≈ 235 years for Neptune-Triton)

### Numerical Implementation
```python
def stability_triangular_point(mu):
    """Analyze stability at triangular Lagrange points"""
    
    # Critical mass ratio
    mu_critical = 0.5 * (1 - np.sqrt(23/27))
    
    if mu > mu_critical:
        print(f"Unstable: μ = {mu:.6f} > μ_critical = {mu_critical:.6f}")
        return None, None, False
    
    # Frequency calculations for stable case
    discriminant = (4 - 27*mu)**2 - 4*(1 - mu)**2
    
    if discriminant < 0:
        print("Error: Negative discriminant in frequency calculation")
        return None, None, False
    
    omega1_sq = 0.5 * (4 - 27*mu + np.sqrt(discriminant))
    omega2_sq = 0.5 * (4 - 27*mu - np.sqrt(discriminant))
    
    omega1 = np.sqrt(omega1_sq)
    omega2 = np.sqrt(omega2_sq)
    
    # Periods in normalized time units
    T1 = 2*np.pi / omega1  # Short period
    T2 = 2*np.pi / omega2  # Long period
    
    return (omega1, T1), (omega2, T2), True
```

## Nonlinear Stability and KAM Theory

### Limitations of Linear Analysis
Linear stability only determines local behavior near equilibrium. For larger amplitudes:
- **Nonlinear resonances**: Can destabilize otherwise stable motion
- **Arnold diffusion**: Slow chaotic transport along resonance chains
- **Secular effects**: Long-term perturbations accumulate

### KAM Theorem Application
For the triangular points when linearly stable:
- **Invariant tori**: Quasi-periodic orbits persist under small perturbations
- **Resonance gaps**: Rational frequency ratios create instability
- **Finite stability**: Stability only guaranteed for finite (but large) amplitudes

### Effective Stability
Even when formally unstable, practical stability can exist:
- **Nekhoroshev estimates**: Exponentially long stability times
- **Lyapunov time**: Characteristic timescale for exponential divergence
- **Practical considerations**: Mission lifetimes vs. instability timescales

## Periodic Orbit Stability

### Floquet Theory
For periodic orbits with period T, stability is determined by the monodromy matrix:
```
M = Φ(T)
```

Where Φ(t) is the fundamental solution matrix of the linearized system.

### Characteristic Multipliers
The eigenvalues of M (characteristic multipliers) determine stability:
- **|λᵢ| < 1**: Stable directions
- **|λᵢ| > 1**: Unstable directions  
- **|λᵢ| = 1**: Marginal stability

### Symmetry Properties
Due to Hamiltonian structure:
- **Reciprocal pairs**: If λ is an eigenvalue, so is 1/λ
- **Complex conjugates**: For real systems
- **Unit determinant**: det(M) = 1

### Stability Index
For planar periodic orbits, the stability index:
```
s = ½|tr(M)|
```

- **s < 1**: Linearly stable (elliptic)
- **s > 1**: Linearly unstable (hyperbolic)
- **s = 1**: Bifurcation boundary (parabolic)

## Lyapunov Stability and Chaos

### Lyapunov Exponents
For long-term stability analysis, compute Lyapunov exponents:
```
λ = lim(t→∞) (1/t) ln(|δr(t)|/|δr(0)|)
```

### Interpretation
- **λ < 0**: Stable (perturbations decay)
- **λ = 0**: Marginal (perturbations neither grow nor decay)
- **λ > 0**: Chaotic (perturbations grow exponentially)

### Numerical Computation
```python
def lyapunov_exponent(initial_state, mu, integration_time, dt=0.01):
    """Compute largest Lyapunov exponent"""
    
    # Initial separation
    delta0 = 1e-12
    state1 = initial_state.copy()
    state2 = initial_state.copy()
    state2[0] += delta0  # Small perturbation in x
    
    times = np.arange(0, integration_time, dt)
    lyap_sum = 0
    count = 0
    
    for i, t in enumerate(times[1:]):
        # Integrate both trajectories
        sol1 = solve_ivp(cr3bp_equations, [times[i], t], state1, 
                        args=(mu,), rtol=1e-12)
        sol2 = solve_ivp(cr3bp_equations, [times[i], t], state2, 
                        args=(mu,), rtol=1e-12)
        
        state1 = sol1.y[:, -1]
        state2 = sol2.y[:, -1]
        
        # Compute separation
        delta = np.linalg.norm(state2 - state1)
        
        # Accumulate logarithmic growth
        lyap_sum += np.log(delta / delta0)
        count += 1
        
        # Renormalize to prevent overflow
        if delta > 1e-6:
            state2 = state1 + delta0 * (state2 - state1) / delta
    
    # Average Lyapunov exponent
    return lyap_sum / (count * dt)
```

## Manifold Theory and Transport

### Stable and Unstable Manifolds
Near hyperbolic points (like collinear Lagrange points):
- **Stable manifold**: Trajectories approaching the equilibrium
- **Unstable manifold**: Trajectories departing from equilibrium

### Heteroclinic and Homoclinic Connections
- **Heteroclinic**: Manifolds connecting different equilibria
- **Homoclinic**: Manifolds returning to the same equilibrium

### Applications
- **Low-energy transfers**: Exploit manifold structures
- **Chaotic mixing**: Understand phase space transport
- **Escape mechanisms**: Pathways for particle loss

## Computational Implementation

### Complete Stability Analysis
```python
def comprehensive_stability_analysis(mu):
    """Complete stability analysis for CR3BP"""
    
    # Find all Lagrange points
    L_points = find_lagrange_points(mu)
    
    results = {}
    
    # Analyze each point
    for name, coords in L_points.items():
        if name in ['L1', 'L2', 'L3']:
            # Collinear point analysis
            eigenvals, omega_z, U_xx, U_yy, U_zz = stability_collinear_point(coords, mu)
            results[name] = {
                'type': 'collinear',
                'stable': False,
                'eigenvalues': eigenvals,
                'out_plane_frequency': omega_z,
                'second_derivatives': (U_xx, U_yy, U_zz)
            }
        else:
            # Triangular point analysis
            freq_data = stability_triangular_point(mu)
            if freq_data[2]:  # If stable
                results[name] = {
                    'type': 'triangular',
                    'stable': True,
                    'short_period_data': freq_data[0],
                    'long_period_data': freq_data[1]
                }
            else:
                results[name] = {
                    'type': 'triangular',
                    'stable': False
                }
    
    return results

def print_stability_summary(results):
    """Print formatted stability analysis results"""
    
    print("Lagrange Point Stability Analysis")
    print("=" * 40)
    
    for name, data in results.items():
        print(f"\n{name} Point:")
        print(f"  Type: {data['type']}")
        print(f"  Stable: {data['stable']}")
        
        if data['type'] == 'collinear':
            print(f"  Out-of-plane frequency: {data['out_plane_frequency']:.6f}")
            print(f"  Eigenvalues: {data['eigenvalues']}")
        elif data['stable'] and data['type'] == 'triangular':
            ω1, T1 = data['short_period_data']
            ω2, T2 = data['long_period_data']
            print(f"  Short-period: ω = {ω1:.6f}, T = {T1:.2f}")
            print(f"  Long-period: ω = {ω2:.6f}, T = {T2:.2f}")
```

## Physical Applications

### Mission Design
- **Halo orbit families**: Periodic orbits around collinear points
- **Station-keeping requirements**: Fuel needed to maintain position
- **Transfer trajectory planning**: Exploit manifold structures

### Natural Dynamics
- **Asteroid stability**: Long-term evolution in multi-body systems
- **Ring particle dynamics**: Stability in planetary ring systems
- **Binary star systems**: Mass transfer and orbital evolution

### Neptune-Triton System
- **Trojan search**: Stability of objects at L₄, L₅
- **Debris evolution**: Long-term fate of impact ejecta
- **Historical dynamics**: Constraints on system evolution

## Advanced Topics

### Parametric Stability
How stability changes with system parameters:
- **Mass ratio variation**: Effects of μ on stability boundaries
- **Perturbation effects**: Impact of additional forces
- **Bifurcation analysis**: Critical parameter values

### High-Order Analysis
Beyond linear approximation:
- **Center manifold theory**: Nonlinear stability near marginal cases
- **Normal forms**: Simplified nonlinear equations
- **Averaging methods**: Long-term behavior analysis

### Numerical Challenges
- **Stiffness**: Wide range of timescales
- **Symplectic integration**: Preserve Hamiltonian structure
- **Error accumulation**: Long-term integration accuracy

## Historical Development

### Classical Period
- **1772**: Lagrange discovers equilibrium points
- **1884**: Lyapunov develops stability theory
- **1892**: Poincaré establishes chaos theory foundations

### Modern Era
- **1954**: Kolmogorov initiates KAM theory
- **1960s**: Computational revolution in stability analysis
- **1990s**: Dynamical systems theory applications to space missions

## References

### Foundational Theory
- Lyapunov, A.M. (1892). "The General Problem of the Stability of Motion"
- Poincaré, H. (1892). *Les méthodes nouvelles de la mécanique céleste*
- Arnold, V.I. (1963). "Proof of a theorem of A.N. Kolmogorov"

### Modern Applications
- Meyer, K.R. & Hall, G.R. (1992). *Introduction to Hamiltonian Dynamical Systems*
- Gómez, G. et al. (2001). *Dynamics and Mission Design Near Libration Points*
- Koon, W.S. et al. (2000). "Low energy transfer to the Moon"

### Computational Methods
- Hairer, E. et al. (2006). *Geometric Numerical Integration*
- Parker, T.S. & Chua, L.O. (1989). *Practical Numerical Algorithms for Chaotic Systems*

---
*Comprehensive framework for understanding stability in three-body dynamics*  
*Essential for mission design and natural system evolution analysis*
