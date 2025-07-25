import numpy as np
from numba import njit, prange
import matplotlib.pyplot as plt

# CRTBP & FTLE Functions
@njit(fastmath=True)
def crtbp_derivatives(x, y, vx, vy, mu):
    x1, y1 = -mu, 0.0
    x2, y2 =  1.0 - mu, 0.0

    r1 = np.hypot(x - x1, y - y1)
    r2 = np.hypot(x - x2, y - y2)

    dvx =  2*vy + x \
         - (1-mu)*(x - x1)/r1**3 \
         -     mu*(x - x2)/r2**3
    dvy = -2*vx + y \
         - (1-mu)*(y - y1)/r1**3 \
         -     mu*(y - y2)/r2**3

    return vx, vy, dvx, dvy

# Runge-Kutta 4th order step for CRTBP
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

    x_new  = x  + (dt/6)*(k1_vx + 2*k2_vx + 2*k3_vx + k4_vx)
    y_new  = y  + (dt/6)*(k1_vy + 2*k2_vy + 2*k3_vy + k4_vy)
    vx_new = vx + (dt/6)*(k1_ax + 2*k2_ax + 2*k3_ax + k4_ax)
    vy_new = vy + (dt/6)*(k1_ay + 2*k2_ay + 2*k3_ay + k4_ay)

    return x_new, y_new, vx_new, vy_new

# Compute FTLE field using Runge-Kutta integration
@njit(parallel=True, fastmath=True)
def compute_ftle(px_vals, py_vals, mu, dt, steps, delta, x0, y0):
    nx, ny = px_vals.size, py_vals.size
    ftle   = np.empty((nx, ny), dtype=np.float64)
    T      = dt * steps

    for i in prange(nx):
        for j in range(ny):
            x1, y1, vx1, vy1 = x0, y0, px_vals[i], py_vals[j]
            x2, y2, vx2, vy2 = x0, y0, px_vals[i] + delta, py_vals[j]
            for _ in range(steps):
                x1, y1, vx1, vy1 = rk4_step(x1, y1, vx1, vy1, mu, dt)
                x2, y2, vx2, vy2 = rk4_step(x2, y2, vx2, vy2, mu, dt)
            d = np.hypot(x2 - x1, y2 - y1)
            ftle[i, j] = np.log(d / delta) / T
    return ftle

# Effective potential for Jacobi constant contours
def effective_potential(x, y, mu):
    r1 = np.hypot(x + mu, y)
    r2 = np.hypot(x - (1-mu), y)
    return 0.5*(x**2 + y**2) + (1-mu)/r1 + mu/r2

def main():
    # System parameters and L4 initial state
    mu = 2.0881498360803843e-4       # Neptune-Triton mass ratio
    x0, y0 = 0.5 - mu, np.sqrt(3)/2  # L4 coordinates

    # Grid and integration settings
    p_max    = 2.0
    n_coarse = 256
    px = np.linspace(-p_max, p_max, n_coarse)
    py = np.linspace(-p_max, p_max, n_coarse)
    dt, steps_T1 = 0.01, 1000
    steps_T2     = 2000
    delta        = 1e-8
    lambda_c     = 0.1

    # Compute FTLE fields
    ftle_T1   = compute_ftle(px, py, mu, dt, steps_T1, delta, x0, y0)
    ftle_T2   = compute_ftle(px, py, mu, dt, steps_T2, delta, x0, y0)
    ftle_diff = ftle_T2 - ftle_T1

    # Adaptive refinement around threshold
    mask = np.abs(ftle_T2 - lambda_c) < 0.02
    if mask.any():
        i_min, i_max = mask.any(axis=1).nonzero()[0][[0, -1]]
        j_min, j_max = mask.any(axis=0).nonzero()[0][[0, -1]]
        px_ref = np.linspace(px[i_min], px[i_max], 512)
        py_ref = np.linspace(py[j_min], py[j_max], 512)
        ftle_ref = compute_ftle(px_ref, py_ref, mu, dt, steps_T2, delta, x0, y0)
    else:
        ftle_ref = None

    # Plot results
    fig, axs = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)

    # FTLE at T2 with C_J contours
    im0 = axs[0].imshow(ftle_T2.T, origin='lower',
                        extent=[-p_max,p_max,-p_max,p_max],
                        cmap='RdYlBu_r', vmin=0, vmax=1.5,
                        interpolation='nearest')
    Omega_0 = effective_potential(x0, y0, mu)
    for CJ in [1.5, 2.0, 2.5]:
        R = np.sqrt(max(0, 2*Omega_0 - CJ))
        axs[0].plot(R*np.cos(np.linspace(0,2*np.pi,400)),
                    R*np.sin(np.linspace(0,2*np.pi,400)), 'k--', lw=1)
    axs[0].set(title='FTLE λ(T₂=2000Δt) with C_J', xlabel='$p_x$', ylabel='$p_y$')
    fig.colorbar(im0, ax=axs[0], label='FTLE')

    # Convergence Delta-Lambda map
    im1 = axs[1].imshow(ftle_diff.T, origin='lower', extent=[-p_max,p_max,-p_max,p_max],
                        cmap='seismic', vmin=-0.5, vmax=0.5, interpolation='nearest')
    axs[1].set(title='FTLE Convergence: Δλ=T₂-T₁', xlabel='$p_x$', ylabel='$p_y$')
    fig.colorbar(im1, ax=axs[1], label='ΔFTLE')

    # Adaptive refinement view
    if ftle_ref is not None:
        im2 = axs[2].imshow(ftle_ref.T, origin='lower', extent=[px_ref[0], px_ref[-1], py_ref[0], py_ref[-1]],
                            cmap='RdYlBu_r', vmin=0, vmax=1.5, interpolation='nearest')
        axs[2].set(title='Adaptive Refinement (λ≈0.1)', xlabel='$p_x$', ylabel='$p_y$')
        fig.colorbar(im2, ax=axs[2], label='FTLE')
    else:
        axs[2].text(0.5, 0.5, 'No refinement region', ha='center', va='center')
        axs[2].axis('off')

    plt.show()

if __name__ == '__main__':
    main()