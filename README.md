# R3BP Neptune-Triton Analysis

## Overview
This repository contains computational analysis of the Neptune-Triton system using the Circular Restricted Three-Body Problem (CR3BP) and related celestial mechanics methods. The project includes Poincaré surface-of-section analysis, stability studies, and perturbation effects including Neptune's J₂ oblateness.

## Project Structure

### `source/CR3BP_Poincare_J2/`
Contains the enhanced CR3BP analysis with optional J₂ perturbation effects:
- **CR3BP_Poincare_j2.py**: Main Poincaré section generator with J₂ perturbation support
- **constants.py**: Physical constants for the Neptune-Triton system
- **parameters.py**: Configuration parameters for simulations
- **plot_PS2_J2.py**: Visualization tools for Poincaré sections

### `source/neptune_triton_models/`
Core Neptune-Triton dynamical models and analysis tools:
- **r3bp/**: Standard three-body problem implementations
- **r3bp_with_perturbations/**: Extended models including various perturbations
- **docs/**: Detailed documentation and analysis results
- **figures/**: Generated plots and visualizations

## Key Features

- **Poincaré Surface-of-Section Analysis**: Generate and visualize phase space structures
- **J₂ Perturbation Effects**: Include Neptune's oblateness in dynamical models
- **Multiple Integration Methods**: Support for both fixed-step and event-detection approaches
- **Stability Analysis**: Lagrange point identification and stability assessment
- **Comparative Studies**: Tools for comparing different modeling approaches

## Getting Started

1. **Install Dependencies**:
   ```bash
   pip install numpy scipy matplotlib numba
   ```

2. **Configure Parameters**:
   Edit the relevant `parameters.py` files to set:
   - Jacobi constant ranges
   - Initial condition sweeps
   - Integration settings
   - Perturbation toggles

3. **Run Simulations**:
   ```bash
   cd source/CR3BP_Poincare_J2/
   python CR3BP_Poincare_j2.py
   ```

4. **Generate Plots**:
   ```bash
   python plot_PS2_J2.py
   ```

## Output

The simulations generate:
- `.dat` files containing Poincaré section crossing data
- Plots showing phase space structure and dynamics
- Analysis results for stability and perturbation effects

Data files are organized by simulation type and parameters for easy comparison and analysis.

## Research Context

This work is part of graduate thesis research investigating the complex dynamics of the Neptune-Triton system. The analysis provides insights into:
- Orbital stability regions
- Chaotic vs. regular motion boundaries  
- Effects of realistic perturbations on idealized models
- Phase space structure and transport mechanisms

## References

Based on methods from:
- Szebehely, V. (1967). *Theory of Orbits*
- Murray, C.D. & Dermott, S.F. (1999). *Solar System Dynamics*
- And contemporary research in celestial mechanics

---
*Author: Blake T. Johnson*  
*Graduate Thesis Project, 2025*
