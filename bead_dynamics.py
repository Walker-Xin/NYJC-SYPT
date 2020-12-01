import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
import time
animation.convert_path: 'C:\Program Files\ImageMagick-7.0.8-Q16\magick.exe'

# Setting constants
omega = np.pi/0.25  # Angular velocity of the loop in radian per second
g = 9.81  # Gravitational accelertaion in meter per second^2
r = 0.1  # Radius of the loop in meter

A = omega**2
B = g / r

phi_0 = 0  # Initial angular displacement from -y-axis anticlockwise of the bead in radian
phi_dot_0 = 0.1  # Initial angular velocity from -y-axis anticlockwise of the bead in radian per second


# Defining differential equations
def dy_dt(y, t):
    # Solving the equation phi_dotdot = A*sin(phi)*cos(phi) - B*sin(phi)
    # Let y be a vector y = [phi_1, phi_2], where phi_1 = phi and phi_2 = phi_dot
    return [y[1], A*np.sin(y[0])*np.cos(y[0]) - B*np.sin(y[0])]


# Generating data with DE
y_0 = [phi_0, phi_dot_0]  # Setting initial values
t_range = np.linspace(0, 10, 1000)
y_range = integrate.odeint(dy_dt, y_0, t_range)
phi_range = y_range[:, 0]
theta_range = -omega * t_range

# Computing gradients
phi_grad = np.diff(phi_range) / np.diff(t_range)
phi_grad = np.append(phi_grad, phi_grad[-1])

# Visualisation
fig, axs = plt.subplots(1, 2, figsize=(12, 7))

axs[0].plot(t_range, phi_range)
axs[0].set_xlabel('t/s')
axs[0].set_ylabel('$\phi$/rad')
axs[1].plot(t_range, phi_grad)
axs[1].set_xlabel('t/s')
axs[1].set_ylabel('$\dot \phi$/rad*s^-1')

plt.text(
    1.01, 0.05, r'$\dot \phi_0$ = {} rad/s'.format(phi_dot_0), transform=plt.gca().transAxes)  # phi_0 text
plt.text(
    1.01, 0.00, r'$\omega$ = {} rad/s'.format(round(omega, 2)), transform=plt.gca().transAxes)  # theta_dot_0 text
plt.text(
    1.01, -0.05, r'$r$ = {} m'.format(round(r, 2)), transform=plt.gca().transAxes)  # Radius text

plt.show()
plt.close()

# 2D animation
fig, ax = plt.subplots(figsize=(8, 8))
ax.axis([-0.25, 0.25, -0.25, 0.25])
ax.set_xlabel('x')
ax.set_ylabel('z')
ax.set_aspect('equal')

x = r * np.sin(phi_range)
y = -r * np.cos(phi_range)

point, = ax.plot([], [], marker='.', color='r')  # Initial bead plot

circle = plt.Circle((0, 0), 0.1, fill=False, lw=0.1)  # Add static loop
ax.add_artist(circle)

time_text = ax.text(
    0.05, 0.95, "Time", transform=ax.transAxes)  # Time text
phi_dot_text = ax.text(
    0.05, 0.90, "phi_dot", transform=ax.transAxes)  # phi_dot text
theta_dot_text = ax.text(
    0.05, 0.85, "$\dot \\theta$ = $\omega$ = {} rad/s".format(round(omega, 2)), transform=ax.transAxes)  # theta_dot text
radius_text = ax.text(
    0.05, 0.80, '$r$ = {} m'.format(round(r, 2)), transform=ax.transAxes)  # Radius text


def animate_2D(i):
    x_coord = x[i]
    y_coord = y[i]
    point.set_data([x_coord], [y_coord])  # Update bead data

    # Update time text
    time_text.set_text('t = {} s'.format(
        round(t_range[i], 2)))  
    # Update phi_dot text
    phi_dot_text.set_text(
        '$\dot \phi$ = {} rad/s'.format(round(phi_grad[i], 2)))
    return point, time_text, phi_dot_text


anim_2D = animation.FuncAnimation(
    fig, animate_2D, frames=range(len(t_range)), interval=50, blit=True)

plt.show()
# Uncomment to save animation
'''start = time.time()
anim_2D.save(r'animation\animation_2D.mp4')
end = time.time()
print('2D saving took {} s'.format(round(end-start, 2)))'''
plt.close()

# 3D animation
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.axes.set_xlim3d(left=-0.15, right=0.15)
ax.axes.set_ylim3d(bottom=-0.15, top=0.15)
ax.axes.set_zlim3d(bottom=-0.15, top=0.15)
ax.grid(False)
# ax.view_init(azim=55, elev=45) # Set viewing angle

# A set of angles for generating the inital loop
circle_angle = np.linspace(0, 2*np.pi, 100)

x = r * np.sin(phi_range) * np.sin(theta_range)
y = r * np.sin(phi_range) * np.cos(theta_range)
z = -r * np.cos(phi_range)

point, = ax.plot([0], [0], [0], marker='.', color='r')  # Initial bead plot

x_circle = np.zeros(np.shape(circle_angle))
y_circle = r * np.cos(circle_angle)
z_circle = r * np.sin(circle_angle)

circle, = ax.plot(
    x_circle, y_circle, z_circle, color='black', lw=0.5)  # Initial loop plot

time_text = ax.text2D(
    0.05, 0.95, "Time", transform=ax.transAxes)  # Time text
phi_dot_text = ax.text2D(
    0.05, 0.90, "phi_dot", transform=ax.transAxes)  # phi_dot text
theta_dot_text = ax.text2D(
    0.05, 0.85, "$\dot \\theta$ = $\omega$ = {} rad/s".format(round(omega, 2)), transform=ax.transAxes)  # theta_dot text
radius_text = ax.text2D(
    0.05, 0.80, '$r$ = {} m'.format(round(r, 2)), transform=ax.transAxes)  # Radius text


def animate_3D(i):
    x_coord = x[i]
    y_coord = y[i]
    z_coord = z[i]
    point.set_data([x_coord], [y_coord])  # Update bead data
    point.set_3d_properties([z_coord], 'z')

    x_circle_new = y_circle * np.sin(theta_range[i])
    y_circle_new = y_circle * np.cos(theta_range[i])
    z_circle_new = z_circle
    circle.set_data(x_circle_new, y_circle_new)  # Update loop data
    circle.set_3d_properties(z_circle_new, 'z')

    # Update time text
    time_text.set_text('t = {} s'.format(
        round(t_range[i], 2)))  
    # Update phi_dot text
    phi_dot_text.set_text(
        '$\dot \phi$ = {} rad/s'.format(round(phi_grad[i], 2)))
    return point, circle, time_text, phi_dot_text


anim_3D = animation.FuncAnimation(fig, animate_3D, frames=range(
    int(len(t_range))), interval=50, blit=True)

plt.show()
# Uncomment to save animation
'''start = time.time()
anim_3D.save(r'animation\animation_3D.mp4')
end = time.time()
print('3D saving took {} s'.format(round(end-start, 2)))'''
plt.close()
