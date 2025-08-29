"""
poincare_trajectory_plot.py
---------------------------------------
Interactive Poincaré Surface-of-Section Explorer for the Neptune–Triton CR3BP (+J₂)

Features:
    • Interactive, browser-based visualization of Poincaré surface-of-section maps.
    • Hover or click any point on the Poincaré map to preview the corresponding trajectory arc in the rotating frame.
    • Identifies and marks the dynamical center (Neptune, Triton, or Lagrange points) about which the selected trajectory is organized.
    • Supports both classical CR3BP and J₂-perturbed models (toggle via parameters.py).
    • Reads pre-generated .dat data files produced by CR3BP_Poincare_J2.py, matching the solver’s naming convention.
    • Automatically adapts to "global" or "highres" mapping modes and perturbed/non-perturbed data folders.
    • All key parameters (Jacobi constant, x₀ sweep, etc.) are set in parameters.py for reproducibility.

Usage:
    1. Set simulation and plotting parameters in parameters.py:
        - C0, DX, XI, XF, mapping_mode, J2_enabled, etc.
    2. Ensure .dat files exist for the chosen parameters/mode in the appropriate data folder.
    3. Run:
        $ python poincare_trajectory_plot.py
    4. The app will:
        - Load all matching .dat files and build a Poincaré map.
        - Open a Dash web app at http://127.0.0.1:8050/
        - Allow you to interactively explore orbits by hovering/clicking on map points.

Author:
    Blake T. Johnson
    Thesis Project
    (c) 2025
"""


from __future__ import annotations
import os, glob, math
import numpy as np

from scipy.integrate import solve_ivp
from scipy.spatial import cKDTree
from scipy.optimize import brentq  # <-- added for robust L1/L2/L3 roots

import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import plotly.graph_objects as go

from parameters import *
from CR3BP_Poincare_j2 import mu  # <-- import mu from your other script

# === Trajectory preview parameters (user adjustable) ===
ORBIT_TSPAN = (-10.0, 10.0)
ORBIT_MAX_STEPS = 5000

# ---- Paths (as you had) ----
project_root = "/Users/blakejohnson/Documents/r3bp_neptune_triton/source/CR3BP_Poincare_J2"
mode_folder = "highres" if str(mapping_mode).lower() == "highres" else "global"
perturb_folder = "perturbed" if J2_enabled else "non_perturbed"

data_folder = os.path.join(project_root, mode_folder, perturb_folder, "data")
results_folder = os.path.join(project_root, mode_folder, perturb_folder, "results")
os.makedirs(results_folder, exist_ok=True)

print(f"[Info] Looking for data in: {data_folder}")
print(f"[Info] Results will be saved in: {results_folder}")

# -------------------- CR3BP utilities --------------------
def r1_r2(x: float, y: float, mu_: float):
    x1, y1 = -mu_, 0.0
    x2, y2 = 1.0 - mu_, 0.0
    dx1, dy1 = x - x1, y - y1
    dx2, dy2 = x - x2, y - y2
    r1 = math.hypot(dx1, dy1)
    r2 = math.hypot(dx2, dy2)
    return r1, r2

# --- Unperturbed (classical) potential ---
def effective_potential_Omega(x: float, y: float, mu_: float) -> float:
    r1, r2 = r1_r2(x, y, mu_)
    return 0.5*(x*x + y*y) + (1.0 - mu_)/r1 + mu_/r2

# --- J2-perturbed potential (planar: z = 0) ---
def effective_potential_Omega_J2(x: float, y: float, mu_: float) -> float:
    try:
        from constants import J2_neptune, R_neptune_meters, a_triton_meters
        J2 = J2_neptune
        R_nd = R_neptune_meters / a_triton_meters
    except Exception:
        J2 = 3.411e-3
        R_nd = 0.0356
    r1, r2 = r1_r2(x, y, mu_)
    Omega = 0.5*(x*x + y*y) + (1.0 - mu_)/r1 + mu_/r2
    Omega += (1.0 - mu_) * J2 * (R_nd**2) / (2.0 * r1**3)
    return Omega

# --- Unperturbed EOM ---
def eom_rotating(t: float, state: np.ndarray, mu_: float) -> np.ndarray:
    x, y, xdot, ydot = state
    r1, r2 = r1_r2(x, y, mu_)
    Ux = x - (1.0 - mu_)*(x + mu_)/r1**3 - mu_*(x - (1.0 - mu_))/r2**3
    Uy = y - (1.0 - mu_)*y/r1**3 - mu_*y/r2**3
    xddot = 2.0*ydot + Ux
    yddot = -2.0*xdot + Uy
    return np.array([xdot, ydot, xddot, yddot])

# --- J2-perturbed EOM (planar) ---
def eom_rotating_J2(t: float, state: np.ndarray, mu_: float) -> np.ndarray:
    try:
        from constants import J2_neptune, R_neptune_meters, a_triton_meters
        J2 = J2_neptune
        R_nd = R_neptune_meters / a_triton_meters
    except Exception:
        J2 = 3.411e-3
        R_nd = 0.0356
    x, y, xdot, ydot = state
    r1, r2 = r1_r2(x, y, mu_)
    Ux = x - (1.0 - mu_)*(x + mu_)/r1**3 - mu_*(x - (1.0 - mu_))/r2**3
    Uy = y - (1.0 - mu_)*y/r1**3 - mu_*y/r2**3
    c = (1.0 - mu_) * 1.5 * J2 * (R_nd**2)
    J2x = -c * (x + mu_) / r1**5
    J2y = -c * y / r1**5
    xddot = 2.0*ydot + (Ux + J2x)
    yddot = -2.0*xdot + (Uy + J2y)
    return np.array([xdot, ydot, xddot, yddot])

# --- Select which to use based on J2_enabled ---
if J2_enabled:
    Omega_func = effective_potential_Omega_J2
    eom_func = eom_rotating_J2
else:
    Omega_func = effective_potential_Omega
    eom_func = eom_rotating

# -------------------- Data loading --------------------
pattern = f"CJ{C0:.5f}_Xi*_DX{DX:.4f}_*_{str(mapping_mode).lower()}_*.dat"
files = sorted(glob.glob(os.path.join(data_folder, pattern)))
if not files:
    print("[Warn] No matching .dat files found.")

FILTER_BY_RANGE = True
if FILTER_BY_RANGE:
    filtered = []
    for f in files:
        try:
            xi_str = [part for part in os.path.basename(f).split("_") if part.startswith("Xi")][0]
            xi_val = float(xi_str[2:])
            if XI <= xi_val <= XF:
                filtered.append(f)
        except Exception:
            continue
    files = filtered

if not files:
    print("[Warn] No files match the specified x-range.")

all_x, all_xdot, index_map = [], [], []
for fi, path in enumerate(files):
    try:
        data = np.loadtxt(path, dtype=float)
    except Exception as e:
        print(f"[Warn] Skipping {path}: {e}")
        continue
    if data.size == 0:
        continue
    if data.ndim == 1:
        data = data[np.newaxis, :]
    nrows = data.shape[0]
    x_vals = data[:, 0].astype(np.float64)
    xdot_vals = data[:, 3].astype(np.float64)
    all_x.append(x_vals)
    all_xdot.append(xdot_vals)
    index_map.extend((fi, ri) for ri in range(nrows))

if len(all_x) > 0:
    X = np.concatenate(all_x)
    Xdot = np.concatenate(all_xdot)
    pts = np.column_stack([X, Xdot])
    tree = cKDTree(pts)
else:
    X = np.array([]); Xdot = np.array([])
    pts = np.empty((0, 2)); tree = None

# -------------------- Trajectory integrator + cache --------------------
class PreviewCache:
    def __init__(self, max_items: int = 128):
        self.max_items = max_items
        self.store: dict[tuple, tuple[np.ndarray, np.ndarray]] = {}
        self.order: list[tuple] = []

    def get(self, key):
        return self.store.get(key, None)

    def put(self, key, value):
        if key in self.store:
            try: self.order.remove(key)
            except ValueError: pass
        elif len(self.order) >= self.max_items:
            oldest = self.order.pop(0)
            self.store.pop(oldest, None)
        self.store[key] = value
        self.order.append(key)

cache = PreviewCache(max_items=128)
cache.store.clear()
cache.order.clear()

def compute_ydot_from_C(x: float, xdot: float, C: float, mu_: float) -> float:
    Omega = Omega_func(x, 0.0, mu_)
    val = 2.0*Omega - C - xdot*xdot
    return math.sqrt(max(val, 0.0))

def integrate_arc(x: float, xdot: float, C: float, mu_: float,
                  t_span=ORBIT_TSPAN, max_steps=ORBIT_MAX_STEPS):
    """Integrate a short arc forward & backward from the y=0 crossing."""
    y = YI
    ydot = compute_ydot_from_C(x, xdot, C, mu_)
    state0 = np.array([x, y, xdot, ydot], dtype=float)

    def fun(t, s):
        return eom_func(t, s, mu_)

    n_eval = max(2, int(max_steps // 2))

    t_eval_fwd = np.linspace(0.0, max(t_span[1], 0.0), n_eval)
    sol_f = solve_ivp(fun, (0.0, t_span[1]), state0, method="DOP853",
                      atol=1e-10, rtol=1e-10, t_eval=t_eval_fwd, dense_output=False)

    t_eval_bwd = np.linspace(0.0, min(t_span[0], 0.0), n_eval)
    sol_b = solve_ivp(fun, (0.0, t_span[0]), state0, method="DOP853",
                      atol=1e-10, rtol=1e-10, t_eval=t_eval_bwd, dense_output=False)

    if not sol_f.success or not sol_b.success:
        raise RuntimeError(f"IVP failed: fwd={sol_f.message} | bwd={sol_b.message}")

    xs = np.concatenate([sol_b.y[0][::-1], sol_f.y[0]]) if sol_b.y.size else sol_f.y[0]
    ys = np.concatenate([sol_b.y[1][::-1], sol_f.y[1]]) if sol_b.y.size else sol_f.y[1]
    if xs.size == 0 or ys.size == 0:
        raise RuntimeError("Integrator returned no samples; increase ORBIT_MAX_STEPS.")
    return xs, ys

# -------------------- L-points + center identification --------------------
def _Ux_total_on_xaxis(x: float, mu_: float) -> float:
    """∂Ω/∂x at y=0 (includes J2 when enabled)."""
    r1 = abs(x + mu_)
    r2 = abs(x - (1.0 - mu_))
    Ux = x - (1.0 - mu_)*(x + mu_)/r1**3 - mu_*(x - (1.0 - mu_))/r2**3
    if J2_enabled:
        try:
            from constants import J2_neptune, R_neptune_meters, a_triton_meters
            J2 = J2_neptune
            R_nd = R_neptune_meters / a_triton_meters
        except Exception:
            J2 = 3.411e-3
            R_nd = 0.0356
        Ux += -(1.0 - mu_) * 1.5 * J2 * (R_nd**2) * (x + mu_) / r1**5
    return Ux

def _bracket_and_root(f, a, b, pieces=400):
    """Find a sign change in [a,b] (subdivide if needed), then brentq."""
    fa = f(a); fb = f(b)
    if fa == 0: return a
    if fb == 0: return b
    if fa*fb < 0:
        return brentq(f, a, b, maxiter=200, xtol=1e-14)
    xs = np.linspace(a, b, pieces+1)
    fs = [f(xi) for xi in xs]
    for i in range(pieces):
        if fs[i] == 0:
            return xs[i]
        if fs[i]*fs[i+1] < 0:
            return brentq(f, xs[i], xs[i+1], maxiter=200, xtol=1e-14)
    # fallback: return the minimizer of |f|
    k = int(np.argmin(np.abs(fs)))
    return float(xs[k])

def compute_colinear_L_points(mu_: float):
    """Return x-coordinates of L1, L2, L3 on y=0 (J2 honored if enabled)."""
    f = lambda x: _Ux_total_on_xaxis(x, mu_)
    x2 = 1.0 - mu_
    x1 = -mu_
    L1 = _bracket_and_root(f, x2 - 0.6, x2 - 1e-8)
    L2 = _bracket_and_root(f, x2 + 1e-8, x2 + 1.5)
    L3 = _bracket_and_root(f, -1.5, x1 - 1e-8)
    return float(L1), float(L2), float(L3)

def compute_triangle_L_points(mu_: float):
    """Classical L4/L5 (J2 shifts ignored for simplicity)."""
    x = 0.5 - mu_
    y = math.sqrt(3)/2.0
    return (x,  y), (x, -y)

def identify_orbit_center(xs: np.ndarray, ys: np.ndarray, mu_: float):
    """Pick the object the arc is 'about' by median distance to refs (robust to outliers)."""
    # Primaries
    refs = {
        "Neptune": (-mu_, 0.0),
        "Triton": (1.0 - mu_, 0.0),
    }
    # L-points
    L1x, L2x, L3x = compute_colinear_L_points(mu_)
    (L4x, L4y), (L5x, L5y) = compute_triangle_L_points(mu_)
    refs.update({
        "L1": (L1x, 0.0), "L2": (L2x, 0.0), "L3": (L3x, 0.0),
        "L4": (L4x, L4y), "L5": (L5x, L5y),
    })

    # Median distance is a good robust score
    best_name, best_xy, best_score = None, None, float("inf")
    for name, (cx, cy) in refs.items():
        d = np.median(np.hypot(xs - cx, ys - cy))
        if d < best_score:
            best_name, best_xy, best_score = name, (cx, cy), d

    return best_name, best_xy, best_score

# -------------------- Dash app --------------------
app = dash.Dash(__name__)

# Left: Poincaré scatter (slightly larger markers for easier hover)
scatter = go.Scattergl(
    x=X, y=Xdot, mode="markers",
    marker=dict(size=4, opacity=0.75),
    hovertemplate="x=%{x:.6f}<br>ẋ=%{y:.6f}<extra></extra>",
)
fig_ps = go.Figure(scatter)
fig_ps.update_layout(
    title=f"Poincaré Map — CJ={C0:.5f}, DX={DX:.4f}, x₀ ∈ [{XI:.3f}, {XF:.3f}]  |  Mode={mode_folder}  |  J₂={'On' if J2_enabled else 'Off'}",
    xaxis_title="x (nondimensional)", yaxis_title="ẋ (nondimensional)",
    hovermode="closest",
)

# Right: orbit preview (empty initially)
fig_orbit = go.Figure()
fig_orbit.update_layout(
    title="Hover or click a point to preview orbit arc",
    xaxis_title="x", yaxis_title="y",
    xaxis=dict(scaleanchor="y", scaleratio=1),
)

app.layout = html.Div([
    html.Div([
        dcc.Graph(id="ps", figure=fig_ps, clear_on_unhover=False),
    ], style={"width": "58%", "display": "inline-block", "verticalAlign": "top"}),
    html.Div([
        dcc.Graph(id="orbit", figure=fig_orbit),
        html.Div(id="status", style={"fontFamily": "monospace", "fontSize": "12px", "marginTop": "8px"})
    ], style={"width": "41%", "display": "inline-block", "marginLeft": "1%"}),
])

# Add click fallback so picking always works
@app.callback(
    Output("orbit", "figure"),
    Output("status", "children"),
    Input("ps", "hoverData"),
    Input("ps", "clickData"),
)
def update_orbit(hoverData, clickData):
    if tree is None or pts.shape[0] == 0:
        return fig_orbit, "No points loaded (KD-tree empty)."

    point = hoverData or clickData
    if point is None or not point.get("points"):
        return fig_orbit, "Hover or click a point to preview orbit arc."

    x_h = float(point["points"][0]["x"])
    xdot_h = float(point["points"][0]["y"])

    # Nearest neighbor in (x, xdot)
    dist, idx = tree.query([x_h, xdot_h], k=1)
    file_idx, row_idx = index_map[int(idx)]

    # Load the exact row
    path = files[file_idx]
    try:
        data = np.loadtxt(path, dtype=float)
        if data.ndim == 1:
            data = data[np.newaxis, :]
    except Exception as e:
        msg = f"Failed to read {os.path.basename(path)}: {e}"
        fig = go.Figure(); fig.update_layout(title=msg)
        return fig, msg

    if not (0 <= row_idx < data.shape[0]):
        msg = f"Row {row_idx} out of bounds for {os.path.basename(path)}."
        fig = go.Figure(); fig.update_layout(title=msg)
        return fig, msg

    row = data[row_idx]
    x_val = float(row[0])
    xdot_val = float(row[3])

    # Cache key now includes preview params to avoid stale arcs
    key = (
        file_idx, row_idx, float(C0),
        float(ORBIT_TSPAN[0]), float(ORBIT_TSPAN[1]),
        int(ORBIT_MAX_STEPS), int(bool(J2_enabled)),
    )
    cached = cache.get(key)
    if cached is None:
        try:
            xs, ys = integrate_arc(x_val, xdot_val, C0, mu,
                                   t_span=ORBIT_TSPAN, max_steps=ORBIT_MAX_STEPS)
            cache.put(key, (xs, ys))
        except Exception as e:
            msg = f"Integration failed: {e}"
            fig = go.Figure(); fig.update_layout(title=msg)
            return fig, msg
    else:
        xs, ys = cached

    # Identify what the orbit is "about" (Neptune, Triton, or L1..L5)
    center_name, (cx, cy), score = identify_orbit_center(xs, ys, mu)

    fig = go.Figure()
    # dot marker + label at the identified center (changed from star to dot, color navy blue)
    fig.add_trace(go.Scatter(
        x=[cx], y=[cy], mode="markers+text",
        marker=dict(size=12, symbol="circle", color="navy"),  # changed symbol and color
        text=[center_name], textposition="bottom center",
        name="center"
    ))
    # the orbit arc (changed color to sky blue)
    fig.add_trace(go.Scatter(
        x=xs, y=ys, mode="lines", name="preview arc",
        line=dict(color="red", width=2)  # changed color and width
    ))
    # the section point (kept as green dot)
    fig.add_trace(go.Scatter(
        x=[x_val], y=[0.0], mode="markers",
        marker=dict(size=8, color="green"), name="section point"
    ))
    fig.update_layout(
        title=f"Orbit preview — file={os.path.basename(path)}  row={row_idx}",
        xaxis_title="x", yaxis_title="y",
        xaxis=dict(scaleanchor="y", scaleratio=1),
        showlegend=False,
    )

    # (Optional) center the view around the chosen center — uncomment to enable:
    # xr = xs.min(), xs.max()
    # yr = ys.min(), ys.max()
    # pad = 0.05 * max(xr[1]-xr[0], yr[1]-yr[0])
    # fig.update_xaxes(range=[min(xr[0], cx)-pad, max(xr[1], cx)+pad])
    # fig.update_yaxes(range=[min(yr[0], cy)-pad, max(yr[1], cy)+pad])

    status = (
        f"Nearest: dist={dist:.2e} | x={x_val:.6f}  ẋ={xdot_val:.6f}  |  CJ={C0:.5f}  |  μ={mu:.6e}  "
        f"| center={center_name} @ ({cx:.6f}, {cy:.6f})  "
        f"| T=({ORBIT_TSPAN[0]:g},{ORBIT_TSPAN[1]:g}) N={ORBIT_MAX_STEPS}"
    )
    return fig, status

print(f"Data folder: {data_folder}")
print(f"Pattern: {pattern}")
print(f"Files found: {files}")

if __name__ == "__main__":
    app.run(debug=True, use_reloader=False)




