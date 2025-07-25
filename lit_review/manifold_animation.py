import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Set up the Lyapunov orbit (simple circle for illustration)
theta = np.linspace(0, 2 * np.pi, 200)
x_orbit = np.cos(theta)
y_orbit = np.sin(theta)

# Define the manifolds as inward (stable) and outward (unstable) spirals
theta_manifold = np.linspace(0, 4 * np.pi, 300)
r_stable = 1 - 0.05 * theta_manifold / (2 * np.pi)
r_unstable = 1 + 0.05 * theta_manifold / (2 * np.pi)

x_stable = r_stable * np.cos(theta_manifold)
y_stable = r_stable * np.sin(theta_manifold)
x_unstable = r_unstable * np.cos(theta_manifold)
y_unstable = r_unstable * np.sin(theta_manifold)

# Set up plot
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.set_aspect('equal')
ax.set_title('Manifold Tubes of a Periodic Orbit')
ax.set_xlabel('x')
ax.set_ylabel('y')

orbit_line, = ax.plot([], [], 'k-', lw=2, label='Periodic Orbit')
stable_line, = ax.plot([], [], 'b--', lw=1.5, label='Stable Manifold')
unstable_line, = ax.plot([], [], 'r--', lw=1.5, label='Unstable Manifold')

ax.legend()

# Animation function
def animate(i):
    n = i + 10
    orbit_line.set_data(x_orbit, y_orbit)
    stable_line.set_data(x_stable[:n], y_stable[:n])
    unstable_line.set_data(x_unstable[:n], y_unstable[:n])
    return orbit_line, stable_line, unstable_line

ani = animation.FuncAnimation(fig, animate, frames=1000, interval=50, blit=True)

plt.show()
