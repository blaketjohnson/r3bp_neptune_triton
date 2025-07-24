import numpy as np
from constants import a_triton_meters, M_neptune, M_triton, G

# Characteristic quantities for the Neptune–Triton system
L = a_triton_meters                      # Distance unit: semi-major axis of Triton
M = M_neptune + M_triton                 # Mass unit: total system mass
mu_nd = M_triton / (M_neptune + M_triton)        # Dimensionless mass parameter
T = np.sqrt(L**3 / (G * M))             # Time unit: from Kepler's third law
V = L / T                                # Velocity unit
E = V**2                                 # Energy per unit mass

# Conversion functions
def to_nondim_length(x_m): return x_m / L
def to_nondim_velocity(v_mps): return v_mps / V
def to_nondim_time(t_s): return t_s / T
def to_nondim_energy(e_j_per_kg): return e_j_per_kg / E

def to_dim_length(x_nd): return x_nd * L
def to_dim_velocity(v_nd): return v_nd * V
def to_dim_time(t_nd): return t_nd * T
def to_dim_energy(e_nd): return e_nd * E

# Print diagnostics
if __name__ == "__main__":
    print("=== Neptune–Triton Nondimensional Scaling ===")
    print(f"Length scale L (a_triton): {L:.3e} m")
    print(f"Mass scale M (Neptune + Triton): {M:.3e} kg")
    print(f"Gravitational constant G: {G:.5e} m^3/kg/s^2")
    print(f"Time scale T: {T:.3e} s")
    print(f"Velocity scale V: {V:.3e} m/s")
    print(f"Energy scale E: {E:.3e} J/kg")
    print(f"Dimensionless mass parameter μ: {mu_nd:.8f}")