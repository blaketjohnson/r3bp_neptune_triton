import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import os
from datetime import datetime

# === Earth-Moon System Constants (Dimensional) ===
M_earth = 5.972e24  # kg
M_moon = 7.342e22   # kg
R_earth = 6378.1e3  # m
R_moon = 1737.4e3   # m
a_earth_moon = 384400e3  # m
T_moon = 27.321661 * 24 * 3600  # s

# === Derived Parameters (Non-Dimensionalized) ===
mu = M_moon / (M_earth + M_moon)
R1_nd = R_earth / a_earth_moon
R2_nd = R_moon / a_earth_moon

# === Control Parameters ===
C = 3.02
XI = 0.18
XF = 0.23
DX = 0.001
tlim_sec = 500000.0
dt_sec = 0.01

def time_to_nd(t_sec):
    return t_sec / T_moon

def gradient_U(x, y, z):
    r1 = np.sqrt((x + mu)**2 + y**2 + z**2)
    r2 = np.sqrt((x - (1 - mu))**2 + y**2 + z**2)
    Ox = x - (1 - mu)*(x + mu)/r1**3 - mu*(x - (1 - mu))/r2**3
    Oy = y - (1 - mu)*y/r1**3 - mu*y/r2**3
    Oz = -(1 - mu)*z/r1**3 - mu*z/r2**3
    return Ox, Oy, Oz

def equations(t, state):
    x, y, z, vx, vy, vz = state
    Ux, Uy, Uz = gradient_U(x, y, z)
    return [vx, vy, vz, Ux + 2 * vy, Uy - 2 * vx, Uz]

def y_cross_event(t, f):
    return f[1]
y_cross_event.direction = 0
y_cross_event.terminal = False

def generate_poincare(C, x0):
    y0, z0 = 0.0, 0.0
    R1 = np.sqrt((x0 + mu)**2)
    R2 = np.sqrt((x0 - (1 - mu))**2)
    arg = -C + x0**2 + 2 * ((1 - mu)/R1 + mu/R2)
    if arg < 0:
        print(f"[SKIP] x0 = {x0:.5f}, imaginary vy0")
        return None
    vy0 = np.sqrt(arg)
    initial_state = [x0, y0, z0, 0.0, vy0, 0.0]

    sol = solve_ivp(
        equations,
        [0, time_to_nd(tlim_sec)],
        initial_state,
        t_eval=np.arange(0, time_to_nd(tlim_sec), time_to_nd(dt_sec)),
        events=[y_cross_event],
        method='DOP853',
        rtol=1e-10,
        atol=1e-12
    )

    crossings = []
    for i in range(1, sol.y.shape[1]):
        if sol.y[1, i-1] * sol.y[1, i] < 0 and sol.y[1, i] > 0:
            xm = 0.5 * (sol.y[:, i] + sol.y[:, i-1])
            crossings.append(xm)

    if not crossings:
        print(f"[NO CROSSINGS] x0 = {x0:.5f}, C = {C:.5f}")
        return None

    data = np.array(crossings)
    fname = f"test_C{C:.5f}_x0{x0:.5f}.dat"
    np.savetxt(fname, data)
    print(f"[DONE] x0 = {x0:.5f}, C = {C:.5f}, points = {len(data)}")
    return fname

def plot_combined_poincare(file_list):
    plt.figure(figsize=(8, 6))
    for fname in file_list:
        data = np.loadtxt(fname)
        if data.ndim == 1:
            data = data[np.newaxis, :]
        plt.scatter(data[:, 0], data[:, 3], s=0.5, label=fname)
    plt.xlabel(r"$x$ (ND)", fontsize=14)
    plt.ylabel(r"$\dot{x}$ (ND)", fontsize=14)
    plt.title(f"Combined Poincar\u00e9 Section at C = {C:.5f}", fontsize=16)
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    print(f"[INFO] mu = {mu:.8f}, R1 = {R1_nd:.5f}, R2 = {R2_nd:.5f}")
    file_list = []
    for x0 in np.arange(XI, XF + DX/2, DX):
        fname = generate_poincare(C, x0)
        if fname:
            file_list.append(fname)

    if file_list:
        plot_combined_poincare(file_list)
    else:
        print("[INFO] No valid crossings to plot.")






