# Import libraries
import numpy as np
from numba import njit, prange
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import plotly.graph_objects as go

# CRTBP & FTLE Functions !(Change to Neptune and Triton)
'''
# Implement the Circular Restricted Three-Body Problem (CRTBP) dynamics and finite-time Lyapunov exponent (FTLE) calculation.

Compute derivatives for CRTBP equations of motion.
Args:
    x, y   - Position coordinates of the test particle
    vx, vy - Velocity components of the test particle
    mu     - Mass ratio of the secondary body

# Returns:
    vx, vy - Velocity components
    ax, ay - Acceleration components due to gravitational forces and Coriolis effects
''' 
@njit(fastmath=True)
def crtbp_derivatives(x, y, vx, vy, mu):
    x1, y1 = -mu, 0.0              # Primary 1 position
    x2, y2 =  1.0 - mu, 0.0        # Primary 2 position
    r1 = np.hypot(x - x1, y - y1)  # Distance to primary 1
    r2 = np.hypot(x - x2, y - y2)  # Distance to primary 2

    # CRTBP equations (include Coriolis and effective potential terms)
    ax =  2*vy + x \
        - (1-mu)*(x - x1)/r1**3 \
        -     mu*(x - x2)/r2**3
    ay = -2*vx + y \
        - (1-mu)*(y - y1)/r1**3 \
        -     mu*(y - y2)/r2**3
    return vx, vy, ax, ay

# RK4 integration step for CRTBP !(change to DOP853)
'''
Perform a single RK4 integration step for CRTBP.
Args:   
    x, y, vx, vy : Current state
    mu           : Mass ratio
    dt           : Time step
Returns:
    x_new, y_new, vx_new, vy_new : Updated state after one RK4 step
'''
@njit(fastmath=True)
def rk4_step(x, y, vx, vy, mu, dt):
    k1_vx, k1_vy, k1_ax, k1_ay = crtbp_derivatives(x, y, vx, vy, mu) 
    x2  = x  + 0.5*dt*k1_vx
    y2  = y  + 0.5*dt*k1_vy
    vx2 = vx + 0.5*dt*k1_ax
    vy2 = vy + 0.5*dt*k1_ay
    k2_vx, k2_vy, k2_ax, k2_ay = crtbp_derivatives(x2, y2, vx2, vy2, mu)
    x3  = x  + 0.5*dt*k2_vx
    y3  = y  + 0.5*dt*k2_vy
    vx3 = vx + 0.5*dt*k2_ax
    vy3 = vy + 0.5*dt*k2_ay
    k3_vx, k3_vy, k3_ax, k3_ay = crtbp_derivatives(x3, y3, vx3, vy3, mu)
    x4  = x  + dt*k3_vx
    y4  = y  + dt*k3_vy
    vx4 = vx + dt*k3_ax
    vy4 = vy + dt*k3_ay
    k4_vx, k4_vy, k4_ax, k4_ay = crtbp_derivatives(x4, y4, vx4, vy4, mu)

    # Weighted sum of increments for RK4
    x_new  = x  + (dt/6)*(k1_vx + 2*k2_vx + 2*k3_vx + k4_vx)
    y_new  = y  + (dt/6)*(k1_vy + 2*k2_vy + 2*k3_vy + k4_vy)
    vx_new = vx + (dt/6)*(k1_ax + 2*k2_ax + 2*k3_ax + k4_ax)
    vy_new = vy + (dt/6)*(k1_ay + 2*k2_ay + 2*k3_ay + k4_ay)
    return x_new, y_new, vx_new, vy_new

# Compute the FTLE field over a grid of initial momenta
'''
Compute the FTLE field over a grid of initial momenta.
Args:
    px_vals, py_vals : Arrays of initial momentum values in x and y directions
    mu              : Mass ratio of the secondary body
    dt              : Time step for integration
    steps           : Number of integration steps
    delta           : Initial separation distance for FTLE calculation
    x0, y0         : Initial position of the test particle

Returns:
    ftle            : 2D array of FTLE values for the momentum grid    
'''
@njit(parallel=True, fastmath=True)
def compute_ftle(px_vals, py_vals, mu, dt, steps, delta, x0, y0):
    nx, ny = px_vals.size, py_vals.size
    ftle  = np.empty((nx, ny), dtype=np.float64)
    T = dt * steps  # Total integration time

    # Loop over momentum grid points
    for i in prange(nx):
        for j in range(ny):
            # Initialize two nearby trajectories in momentum space
            x1, y1, vx1, vy1 = x0, y0, px_vals[i], py_vals[j]
            x2, y2, vx2, vy2 = x0, y0, px_vals[i]+delta, py_vals[j]

            # Integrate both trajectories forward in time
            for _ in range(steps):
                x1, y1, vx1, vy1 = rk4_step(x1, y1, vx1, vy1, mu, dt)
                x2, y2, vx2, vy2 = rk4_step(x2, y2, vx2, vy2, mu, dt)

            # Compute separation after integration
            d = np.hypot(x2 - x1, y2 - y1)
            # FTLE formula: normalized log growth rate of separation
            ftle[i, j] = np.log(d/delta)/T

    return ftle

# Integrate a single orbit trajectory for visualization
'''
Integrate a single orbit trajectory for visualization.
Args:
    x0, y0, vx0, vy0 : Initial state (position and velocity)
    mu               : Mass ratio
    dt               : Time step for integration
    steps            : Number of integration steps

Returns:
    xs, ys          : Arrays of x and y positions over the trajectory
'''
def integrate_orbit(x0, y0, vx0, vy0, mu, dt, steps):
    xs = np.empty(steps+1)
    ys = np.empty(steps+1)
    x, y, vx, vy = x0, y0, vx0, vy0
    xs[0], ys[0] = x, y
    for k in range(1, steps+1):
        x, y, vx, vy = rk4_step(x, y, vx, vy, mu, dt)
        xs[k], ys[k] = x, y
    return xs, ys

# FTLE Field Precomputation
# Set up parameters for the CRTBP and FTLE calculation
mu = 2.0881498360803843e-4                 # Mass ratio !(Earth-Moon system, change to Neptune and Triton)
x0, y0 = 0.5 - mu, np.sqrt(3)/2            # Initial position (L4 point)
p_max, n_coarse = 2.0, 1024                # Range and resolution for momentum grid
px = np.linspace(-p_max, p_max, n_coarse)  # Momentum in x direction grid
py = np.linspace(-p_max, p_max, n_coarse)  # Momentum in y direction grid
dt, steps_T2, delta = 0.01, 10000, 1e-8    # Integration time step, steps, FTLE delta

print("Computing FTLE field…")

# Compute FTLE field over the momentum grid
ftle = compute_ftle(px, py, mu, dt, steps_T2, delta, x0, y0).T
print("Done.")

# Create Dash app for interactive FTLE and orbit visualization
app = dash.Dash(__name__)

# Create FTLE heatmap figure using Plotly
fig_heat = go.Figure(go.Heatmap(z=ftle, x=px, y=py,
                                colorscale='RdYlBu_r', zmin=0, zmax=1.5, colorbar=dict(title='FTLE')))
fig_heat.update_layout(title="FTLE λ(T₂=2000 Δt)", hovermode='closest')

# Create empty orbit figure for updates
fig_orbit = go.Figure(go.Scatter(x=[], y=[], mode='lines'))
fig_orbit.update_layout(title="Hover over heatmap to see orbit", xaxis_title='x', yaxis_title='y')

# FTLE heatmap and orbit plot
app.layout = html.Div([html.Div(dcc.Graph(id='heatmap', figure=fig_heat), style={'width':'60%', 'display':'inline-block'}),
                       html.Div(dcc.Graph(id='orbit',   figure=fig_orbit), style={'width':'38%', 'display':'inline-block', 'vertical-align':'top'})])

# Update orbit plot when user hovers over FTLE heatmap
@app.callback(Output('orbit', 'figure'), Input('heatmap', 'hoverData'))
def update_orbit(hoverData):
    if hoverData is None:
        return fig_orbit
    # Find nearest grid point to hovered location
    x_h = hoverData['points'][0]['x']
    y_h = hoverData['points'][0]['y']
    i = np.abs(px - x_h).argmin()
    j = np.abs(py - y_h).argmin()
    xs, ys = integrate_orbit(x0, y0, px[i], py[j], mu, dt, steps_T2)
    fig = go.Figure(go.Scatter(x=xs, y=ys, mode='lines'))
    fig.update_layout(title=f"Orbit at (pₓ={px[i]:.3f}, pᵧ={py[j]:.3f})",
                      xaxis_title='x', yaxis_title='y')
    return fig

# Run the Dash app for interactive visualization
if __name__ == '__main__':
    app.run(debug=True)

