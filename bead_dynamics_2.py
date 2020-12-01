import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3

# Setting constants
omega = np.pi/0.25  # Angular velocity of the loop in radian per second
g = 9.81  # Gravitational accelertaion in meter per second^2
r = 0.1  # Radius of the loop in meter

A = omega**2
B = g / r

phi_0 = 0  # Initial angular displacement from -y-axis anticlockwise of the bead in radian
phi_dot_0 = 20  # Initial angular velocity from -y-axis anticlockwise of the bead in radian per second


# Defining differential equations
def dy_dt(y, t):
    # Solving the equation phi_dotdot = A*sin(phi)*cos(phi) - B*sin(phi)
    # Let y be a vector y = [phi_1, phi_2], where phi_1 = phi and phi_2 = phi_dot
    return [y[1], A*np.sin(y[0])*np.cos(y[0]) - B*np.sin(y[0])]


# Generating data
y_0 = [phi_0, phi_dot_0]  # Setting initial values
t_range = np.linspace(0, 20, 1000)
y_range = integrate.odeint(dy_dt, y_0, t_range)
phi_range = y_range[:, 0]

# Visualisation
plt.xlabel('t/s')
plt.ylabel('$\phi$/rad')

plt.text(0.85, 0.15, r'$\dot \phi_0$ = {} rad'.format(phi_dot_0),
         {'fontsize': 10}, transform=plt.gca().transAxes)
plt.text(0.85, 0.10, r'$\omega$ = {} rad/s'.format(round(omega, 2)),
         {'fontsize': 10}, transform=plt.gca().transAxes)
plt.text(0.85, 0.05, r'$r$ = {} m'.format(round(r, 2)), {
         'fontsize': 10}, transform=plt.gca().transAxes)

plt.plot(t_range, phi_range)
plt.show()
plt.close()

# 2D animation
fig, ax = plt.subplots(figsize=(8,8))
ax.axis([-0.25,0.25,-0.25,0.25])
ax.set_xlabel('x')
ax.set_ylabel('z')
ax.set_aspect('equal')

point, = ax.plot([], [], marker='.', color='r')

circle = plt.Circle((0,0), 0.1, fill=False, lw=0.1)
ax.add_artist(circle)

x = r * np.sin(phi_range)

y = -r * np.cos(phi_range)

time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def animate_2D(i):
    x_coord = x[i]
    y_coord = y[i]
    point.set_data([x_coord], [y_coord])

    time_text.set_text('t = {} s'.format(round(t_range[i], 2)))
    return point, time_text


anim = animation.FuncAnimation(fig, animate_2D, frames=range(len(x)), interval=10, blit=True)

plt.show()
plt.close()

# 3D animation
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.axes.set_xlim3d(left=-0.15, right=0.15)
ax.axes.set_ylim3d(bottom=-0.15, top=0.15)
ax.axes.set_zlim3d(bottom=-0.15, top=0.15)
ax.grid(False)
# ax.view_init(azim=55, elev=45)

point, = ax.plot([0], [0], [0], marker='.', color='r')

theta_range = -omega * t_range

x = r * np.sin(phi_range) * np.sin(theta_range)

y = r * np.sin(phi_range) * np.cos(theta_range)

z = -r * np.cos(phi_range)

x_circle = np.zeros(np.shape(theta_range))

y_circle = r * np.cos(theta_range)

z_circle = r * np.sin(theta_range)

circle, = ax.plot(x_circle, y_circle, z_circle, color='black', lw=0.1)

time_text = ax.text2D(0.05, 0.95, "2D Text", transform=ax.transAxes)


def animate_3D(i):
    x_coord = x[i]
    y_coord = y[i]
    z_coord = z[i]
    point.set_data([x_coord], [y_coord])
    point.set_3d_properties([z_coord], 'z')
    x_circle_new = y_circle * np.sin(theta_range[i])
    y_circle_new = y_circle * np.cos(theta_range[i])
    z_circle_new = z_circle
    circle.set_data(x_circle_new, y_circle_new)
    circle.set_3d_properties(z_circle_new, 'z')
    time_text.set_text('t = {} s'.format(round(t_range[i], 2)))
    return point, circle, time_text


anim = animation.FuncAnimation(fig, animate_3D, frames=range(len(t_range)), interval=1, blit=True)

plt.show()

writergif = animation.PillowWriter(fps=30) 

anim.save('ok.gif', writer=writergif)

plt.close()