import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Angle parameter
theta = np.linspace(0, 6 * np.pi, 600)

# Smooth transition in x from orbit A to orbit B
x_shift = np.linspace(-1.5, 1.5, len(theta))

# Radius function: spiral expanding then contracting
r_base = 0.5 + 0.03 * np.sin(0.5 * theta) + 0.02 * theta

# Generate smooth heteroclinic spiral curve
x_path = x_shift + r_base * np.cos(theta)
y_path = r_base * np.sin(theta)

# Orbit centers
center_A = (-1.5, 0)
center_B = (1.5, 0)

# Setup plot
fig, ax = plt.subplots(figsize=(7, 4))
ax.set_xlim(-3.5, 3.5)
ax.set_ylim(-2, 2)
ax.set_aspect('equal')
ax.set_title("Smooth Heteroclinic Connection")
ax.set_xlabel("x")
ax.set_ylabel("y")

# Draw the periodic orbits
circle_A = plt.Circle(center_A, 0.5, color='k', fill=False, lw=2, label="Orbit A (L₁)")
circle_B = plt.Circle(center_B, 0.5, color='k', fill=False, lw=2, label="Orbit B (L₂)")
ax.add_patch(circle_A)
ax.add_patch(circle_B)

# Animated line and moving point
line, = ax.plot([], [], 'r-', lw=2, label="Heteroclinic Trajectory")
point, = ax.plot([], [], 'ro', markersize=5)
ax.legend(loc='upper center')

# Animation function
def animate(i):
    line.set_data(x_path[:i], y_path[:i])
    if i > 0:
        point.set_data([x_path[i - 1]], [y_path[i - 1]])
    else:
        point.set_data([], [])
    return line, point

ani = animation.FuncAnimation(fig, animate, frames=len(x_path), interval=20, blit=True)

plt.show()

