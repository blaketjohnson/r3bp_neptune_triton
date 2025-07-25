import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Parameters
theta = np.linspace(0, 6 * np.pi, 600)  # angle
r_out = 1 + 0.05 * theta                # spiraling out (unstable)
r_in = r_out[::-1]                      # spiraling back in (stable)

# Combine into one homoclinic path (out + in)
x_out = r_out * np.cos(theta)
y_out = r_out * np.sin(theta)
x_in = r_in * np.cos(theta)
y_in = r_in * np.sin(theta)

x_full = np.concatenate([x_out, x_in])
y_full = np.concatenate([y_out, y_in])

# Setup figure
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_xlim(-3, 3)
ax.set_ylim(-3, 3)
ax.set_aspect('equal')
ax.set_title("Homoclinic Orbit Animation")
ax.set_xlabel("x")
ax.set_ylabel("y")

# Lyapunov orbit (reference periodic orbit)
circle = plt.Circle((0, 0), 1.0, color='k', fill=False, linewidth=2, label="Periodic Orbit")
ax.add_patch(circle)

# Animated trajectory
line, = ax.plot([], [], 'r-', lw=2, label="Homoclinic Orbit")
point, = ax.plot([], [], 'ro', markersize=5)

ax.legend()

# Animation function
def animate(i):
    line.set_data(x_full[:i], y_full[:i])
    if i > 0:
        point.set_data([x_full[i-1]], [y_full[i-1]])
    else:
        point.set_data([], [])
    return line, point



ani = animation.FuncAnimation(fig, animate, frames=len(x_full), interval=20, blit=True)

plt.show()
