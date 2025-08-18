# Neptune–Triton Restricted Three-Body Problem (CR3BP + J₂)

[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Numerical experiments in the **Circular Restricted Three-Body Problem (CR3BP)** with applications to the Neptune–Triton system.  
Implements both **classical CR3BP** and an extended model with **Neptune’s J₂ oblateness perturbation**, producing Poincaré maps, stability regions, and trajectories.

---

## ✨ Features
- Classical **CR3BP equations of motion**
- **J₂-perturbed CR3BP** for Neptune’s oblateness
- Automated **Jacobi constant calculation**
- Root solvers for **Lagrange points (L1–L3)**
- Batch generation of **Poincaré surfaces of section**
- Parallelized integration (DOP853 with event detection)

---

## 📊 Example Results

### Poincaré Surfaces of Section
- Classical CR3BP (stability islands appear as expected)  
<p align="center">
  <img src="source/CR3BP_Poincare_J2/highres/non_perturbed/results/Poincare_C3.01000_DX0.1000_highres.png" alt="Poincaré surface (CR3BP)" width="500"/>
</p>

- CR3BP + J₂ perturbation (oblateness slightly modifies stability structures)  
<p align="center">
  <img src="source/CR3BP_Poincare_J2/highres/non_perturbed/results/Poincare_C3.01430_DX0.0005_highres.png" alt="Poincaré surface (CR3BP+J₂)" width="500"/>
</p>

---

## 🚀 Quick Start

Clone the repo and install dependencies:
```bash
git clone git@github.com:blaketjohnson/r3bp_neptune_triton.git
cd r3bp_neptune_triton
pip install -r requirements.txt
```

Run a sample Poincaré map:
```bash
python scripts/run_poincare.py
```

Run Jacobi constant solver:
```bash
python scripts/jacobi_solver.py
```

---

## 🛠️ Dependencies
- Python 3.10+
- NumPy
- SciPy
- Matplotlib
- PyTest (optional, for unit tests)

Install all requirements with:
```bash
pip install -r requirements.txt
```

---

## 📂 Repository Structure
```
r3bp_neptune_triton/
├── src/
│   ├── equations/       # CR3BP and J₂ equations of motion
│   ├── rootfinding/     # L1–L3 solvers
│   ├── poincare/        # Poincaré event handlers + plotting
│   ├── utils/           # constants, scaling, wrappers
├── scripts/
│   ├── run_poincare.py
│   ├── jacobi_solver.py
│   └── plot_utils.py
├── docs/
│   ├── poincare_cr3bp.png
│   ├── poincare_j2.png
└── tests/
    ├── test_equations.py
    └── test_poincare.py
```

---

## 📚 Background
This project builds on classical CR3BP methods and extends them with Neptune’s oblateness.  
It demonstrates numerical analysis of stability regions and surface of section methods for dynamical astronomy.  
These techniques are also relevant to **mission trajectory design, orbital stability analysis, and spacecraft dynamics** in planetary systems.

---

## 📜 License
MIT License. See LICENSE for details.
