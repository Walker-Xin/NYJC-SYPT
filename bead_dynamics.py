import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
import time

import warnings
# warnings.simplefilter('error') # Uncomment to make script stop should a warning is raised

# Setting parameters
m = 0.2  # Mass of the bead in kilogram
m_0 = 0.5  # Mass of the loop in kilogram
r = 0.1  # Radius of the loop in meter
phi_0 = 30/180 * np.pi # Initial angular displacement from +y-axis clockwise of the bead in radian
omega = np.pi/1  # Initial angular velocity of the loop in radian per second
g = 9.81  # Gravitational accelertaion in meter per second^2

I = 0.5 * m_0 * (r**2)  # Moment of inertia of the loop about the diameter
M = (m * (np.sin(phi_0)**2) * (r**2) + I) * omega  # Initial angular momentum of the system
E = 0.5 * m * (np.sin(phi_0)**2) * (r**2) * (omega**2) + 0.5 * I * (omega **2) + m * g * r * (1 + np.cos(phi_0))  # Initial energy of the system

A = g / r
B = M**2
C = m * (r**2)
# Defining integral functions


def t_integral(phi):
    # Expressing time as a function of phi
    c_1 = 2 * E / (m*(r**2))
    c_2 = M**2
    c_3 = (m**2) * (r**4)
    c_4 = I * m * (r**2)
    c_5 = 2 * g / r
    integral = 1 / np.sqrt(c_1 - c_2 / (c_3 * (np.sin(phi)**2) + c_4) - c_5 * (1 + np.cos(phi)))
    return integral


def theta_integral(phi):
    # Expressing theta as a function of phi
    c_1 = 2 * E / (m*r*r)
    c_2 = M**2
    c_3 = (m**2) * (r**4)
    c_4 = I * m * (r**2)
    c_5 = 2 * g / r
    integral_1 = M / ((m * (r**2) * (np.sin(phi)**2)) + I)
    integral_2 = 1 / np.sqrt(c_1 - c_2 / (c_3 * (np.sin(phi)**2) + c_4) - c_5 * (1 + np.cos(phi)))
    integral = integral_1 * integral_2
    return integral


# Defining differential equations
def dy_dt(y, t):
    # Solving the equation phi_dotdot = A*sin(phi) - B*sin(phi)*cos(phi)/(C*sin^2(phi) + I)^2
    # Let y be a vector y = [phi_1, phi_2], where phi_1 = phi and phi_2 = phi_dot
    return [y[1], A*np.sin(y[0]) - B*np.sin(y[0])*np.cos(y[0])/((C*np.sin(y[0])*np.sin(y[0]) + I)**2)]


# Computing DE
y_0 = [phi_0, 0]
t_range = np.linspace(0, 10, 1000)
y_range = integrate.odeint(dy_dt, y_0, t_range)
phi_range = y_range[:, 0]

# Computing integrals
theta_range = np.zeros(phi_range.size)
for i in range(phi_range.size):
    theta_range[i], error = integrate.quad(theta_integral, phi_0, phi_range[i])

# Computing gradients
theta_grad = np.diff(theta_range) / np.diff(t_range)
phi_grad = np.diff(phi_range) / np.diff(t_range)
theta_grad = np.append(theta_grad, theta_grad[-1])
phi_grad = np.append(phi_grad, phi_grad[-1])

# Visualisation
fig, axs = plt.subplots(2, 2, figsize=(12, 7))

axs[0][0].plot(t_range, theta_range)
axs[0][0].set_yticks(np.arange(np.pi, 3*np.pi, np.pi))
axs[0][0].set_xlabel('t/s')
axs[0][0].set_ylabel('$\\theta$/rad')
axs[0][1].plot(t_range, phi_range)
axs[0][1].set_yticks(np.arange(np.pi, 3*np.pi, np.pi))
axs[0][1].set_xlabel('t/s')
axs[0][1].set_ylabel('$\phi$/rad')
axs[1][0].plot(t_range, theta_grad)
axs[1][0].set_xlabel('t/s')
axs[1][0].set_ylabel('$\dot \\theta$/rad*s^-1')
axs[1][1].plot(t_range, phi_grad)
axs[1][1].set_xlabel('t/s')
axs[1][1].set_ylabel('$\dot \phi$/rad*s^-1')

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
# ax.view_init(azim=0, elev=90)

point, = ax.plot([0], [0], [0], marker='.', color='r')

theta_range = -omega * t_range

x = r * np.sin(phi_range) * np.sin(theta_range)

y = r * np.sin(phi_range) * np.cos(theta_range)

z = r * np.cos(phi_range)

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
    time_text.set_text('t = {} s'.format(round(t_range[i], 2))) # Update time text
    return point, circle, time_text


anim_3D = animation.FuncAnimation(fig, animate_3D, frames=range(len(x)), interval=50, blit=True)

plt.show()
start = time.time()
anim_3D.save('animation_3D_free_loop.mp4')
end = time.time()
print('3D Saving took {} s'.format(round(end-start, 2)))
plt.close()