# CR3BP with J₂ Perturbation Analysis

## Overview
This module implements the Circular Restricted Three-Body Problem (CR3BP) for the Neptune-Triton system with optional inclusion of Neptune's J₂ oblateness perturbation. It generates Poincaré surface-of-section data to analyze phase space structure and orbital dynamics.

## Key Files

### Core Scripts
- **`CR3BP_Poincare_j2.py`**: Main simulation engine
  - Generates Poincaré section crossings (y=0, ẏ>0)
  - Supports both "global" (fixed-step) and "highres" (event-detection) modes
  - Optional J₂ perturbation effects
  - Automatic trajectory termination for escapes and collisions

- **`plot_PS2_J2.py`**: Visualization and analysis tool
  - Creates publication-quality Poincaré section plots
  - Supports multiple data file overlays
  - Customizable plotting parameters

### Configuration Files
- **`constants.py`**: Physical constants for Neptune-Triton system
  - Mass ratios, orbital parameters
  - J₂ coefficient and Neptune radius
  - Unit conversions and normalizations

- **`parameters.py`**: Simulation configuration
  - Jacobi constant ranges and resolution
  - Initial condition sweep parameters  
  - Integration settings and tolerances
  - Safety distances and termination criteria
  - J₂ perturbation toggle

## Usage

### 1. Configure Simulation Parameters
Edit `parameters.py` to set:
```python
# Jacobi constant range
c_start = 3.0000
c_end = 3.0200
c_step = 0.0002

# Initial x-coordinate sweep  
x0_start = 0.1
x0_end = 1.0
x0_step = 0.1

# Mapping mode
mapping_mode = "highres"  # or "global"

# J₂ perturbation
include_j2 = True
```

### 2. Run Simulation
```bash
python CR3BP_Poincare_j2.py
```

### 3. Generate Plots
```bash
python plot_PS2_J2.py
```

## Output Structure

### Data Files
Poincaré section data is saved as `.dat` files in mode-specific directories:
- `Poincare_data_global/` - Fixed-step integration results
- `Poincare_data_highres/` - Event-detection results

Filename format: `PY-C{jacobi_constant}Xi{x0_value}_{mode}.dat`

Example: `PY-C3.01000Xi0.50000_highres.dat`

### Data Format
Each `.dat` file contains space-separated columns:
```
x_crossing  vx_crossing  time_crossing
```

### Plots
Generated plots show:
- Poincaré section structure in (x, vₓ) coordinates
- Different colors/symbols for different Jacobi constants
- Regular vs. chaotic motion regions
- Perturbation effects when J₂ is enabled

## Integration Modes

### Global Mode (`mapping_mode = "global"`)
- Fixed time-step integration with midpoint interpolation
- Faster computation, lower precision crossings
- Good for exploratory analysis and broad parameter sweeps

### High-Resolution Mode (`mapping_mode = "highres"`)
- Event-detection for exact y=0 crossings
- Higher precision, slower computation
- Best for detailed phase space analysis

## J₂ Perturbation Effects

When `include_j2 = True`, the equations of motion include Neptune's oblateness:
- Modifies effective potential and force field
- Creates additional resonance structures
- Affects stability of periodic orbits
- Most significant for orbits close to Neptune

## Termination Criteria

Trajectories are automatically terminated when:
1. **Escape**: Distance exceeds Hill radius (≈ 4.7 Neptune-Triton distance)
2. **Neptune Collision**: Distance to Neptune < `min_distance_neptune`
3. **Triton Collision**: Distance to Triton < `min_distance_triton`
4. **Maximum Time**: Integration time exceeds `t_max`

## Performance Tips

- Use "global" mode for parameter space exploration
- Use "highres" mode for final publication plots
- Adjust `c_step` and `x0_step` based on desired resolution vs. runtime
- Enable J₂ perturbation selectively (increases computation time ~20%)

## Physical Interpretation

The Poincaré sections reveal:
- **Elliptic islands**: Quasi-periodic, stable motion
- **Chaotic seas**: Irregular, unpredictable trajectories  
- **Separatrices**: Boundaries between different motion types
- **Resonance chains**: Networks of periodic orbit families

J₂ perturbations typically:
- Shift resonance locations
- Create new periodic orbit families
- Modify island-chain structures
- Affect long-term stability

---
*Part of Neptune-Triton thesis research project*  
*Author: Blake T. Johnson, 2025*
