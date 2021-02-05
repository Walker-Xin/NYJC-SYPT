import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
import time
import sys
import csv

m = 0.1
g = 9.81
r = 0.05

D_psi_0 = -0.3*np.pi
D_theta_0 = 0
D_phi_0 = 0.3*np.pi
psi_0 = 0
theta_0 = np.radians(80)
phi_0 = 0
x_0 = 0
y_0 = 0
z_0 = 0

D_phi_0 = -((3/2)*np.sin(theta_0)*(D_psi_0**2) +g/r)/(2*np.tan(theta_0)*(D_psi_0))

t_max=10 # Simulation time in seconds
iterations=1000 # Total number of iterations
t_step=t_max/iterations # Simulation time step
print('Producing simulation with {}s between frames for {} seconds...'.format(t_step, t_max))
t_range = np.linspace(0, t_max, iterations)

w_0 = [D_psi_0, D_theta_0, D_phi_0, psi_0, theta_0, phi_0, x_0, y_0]


def dy_dt_matrix(w, t):
    D_psi, D_theta, D_phi, psi, theta, phi, x, y = w
    
    M = np.array([[0.5*np.sin(theta), 0,   0, 0, 0, 0, 0, 0],
                 [0,                  1.5, 0, 0, 0, 0, 0, 0],
                 [2*np.cos(theta),    0,   2, 0, 0, 0, 0, 0],
                 [0,                  0,   0, 1, 0, 0, 0, 0],
                 [0,                  0,   0, 0, 1, 0, 0, 0],
                 [0,                  0,   0, 0, 0, 1, 0, 0],
                 [0,                  0,   0, 0, 0, 0, 1, 0],
                 [0,                  0,   0, 0, 0, 0, 0, 1]])
    
    f = np.array([[D_theta*D_phi*np.cos(theta) + D_theta*D_phi], 
                  [-(3/2)*(D_psi**2)*np.sin(theta)*np.cos(theta) - 2*D_phi*D_psi*np.sin(theta) - (g/r)*np.cos(theta)],
                  [3*D_theta*D_psi*np.sin(theta)],
                  [D_theta],
                  [D_phi],
                  [D_psi],
                  [-r*(D_psi*np.cos(psi)*np.cos(theta) - D_theta*np.sin(psi)*np.sin(theta) + D_phi*np.cos(psi))],
                  [-r*(D_psi*np.sin(psi)*np.cos(theta) + D_theta*np.cos(psi)*np.sin(theta) + D_phi*np.sin(psi))]])
    
    dy_dt = np.linalg.solve(M, f)
    dy_dt = np.squeeze(dy_dt)
    
    return dy_dt


def dy_dt(w, t):
    D_psi, D_theta, D_phi, psi, theta, phi, x, y = w

    DD_psi = (D_theta*D_phi*np.cos(theta) + D_theta*D_phi)/(0.5*np.sin(theta))
    DD_theta = (-(3/2)*(D_psi**2)*np.sin(theta)*np.cos(theta) - 2*D_phi*D_psi*np.sin(theta) - (g/r)*np.cos(theta))/(1.5)
    DD_phi = (3*D_theta*D_psi*np.sin(theta) - 2*np.cos(theta)*DD_psi)/(2)
    D_x = -r*(D_psi*np.cos(psi)*np.cos(theta) - D_theta*np.sin(psi)*np.sin(theta) + D_phi*np.cos(psi))
    D_y = -r*(D_psi*np.sin(psi)*np.cos(theta) + D_theta*np.cos(psi)*np.sin(theta) + D_phi*np.sin(psi))

    dy_dt = [DD_psi, DD_theta, DD_phi, D_psi, D_theta, D_phi, D_x, D_y]

    return dy_dt


w_range = integrate.odeint(dy_dt, w_0, t_range)
psi_dot_range = w_range[:, 0]
theta_dot_range = w_range[:, 1]
phi_dot_range = w_range[:, 2]
psi_range = w_range[:, 3]
theta_range = w_range[:, 4]
phi_range = w_range[:, 5]
x_range = w_range[:, 6]
y_range = w_range[:, 7]
z_range = np.zeros(np.size(x_range))

fig, axs = plt.subplots(2, 3, figsize=(12, 7))

color = 'tab:red'
axs[0][0].plot(t_range, theta_range, color=color)
axs[0][0].set_xlabel('t/s')
axs[0][0].set_ylabel('$\\theta$/rad')
color = 'coral'
axs[1][0].plot(t_range, theta_dot_range, color=color)
axs[1][0].set_xlabel('t/s')
axs[1][0].set_ylabel('$\dot \\theta$/rad$\cdot$s$^-1$')

color = 'tab:red'
axs[0][1].plot(t_range, phi_range, color=color)
axs[0][1].set_xlabel('t/s')
axs[0][1].set_ylabel('$\\phi$/rad')
color = 'coral'
axs[1][1].plot(t_range, phi_dot_range, color=color)
axs[1][1].set_xlabel('t/s')
axs[1][1].set_ylabel('$\dot \\phi$/rad$\cdot$s$^-1$')

color = 'tab:red'
axs[0][2].plot(t_range, psi_range, color=color)
axs[0][2].set_xlabel('t/s')
axs[0][2].set_ylabel('$\\psi$/rad')
color = 'coral'
axs[1][2].plot(t_range, psi_dot_range, color=color)
axs[1][2].set_xlabel('t/s')
axs[1][2].set_ylabel('$\dot \\psi$/rad$\cdot$s$^-1$')


plt.text(
    1.01, 0.075, '$m$ = {} kg'.format(m, 2), transform=plt.gca().transAxes)  # m text
plt.text(
    1.01, 0.00, '$\\theta_0$ = {} rad'.format(round(theta_0, 2)), transform=plt.gca().transAxes)  # theta_0 text
plt.text(
    1.01, -0.075, '$\\phi_0$ = {} rad'.format(round(phi_0, 2)), transform=plt.gca().transAxes)  # phi_0 text
plt.text(
    1.01, -0.15, '$\\psi_0$ = {} rad'.format(round(psi_0, 2)), transform=plt.gca().transAxes)  # psi_0 text

plt.show()
plt.close()

fig, axs = plt.subplots(1, 1, figsize=(12, 12))

color = 'tab:red'
axs.plot(x_range, y_range, color=color)

plt.show()
plt.close()

# 3D animation
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.axes.set_xlim3d(left=2*np.min(x_range)-r, right=2*np.max(x_range)+r)
ax.axes.set_ylim3d(bottom=2*np.min(y_range)-r, top=2*np.max(y_range)+r)
ax.axes.set_zlim3d(bottom=-10*(2*r), top=10*(2*r))
ax.set_box_aspect((1, 1, 1), zoom=1)
ax.view_init(azim=0, elev=0) # Set viewing angle

# A set of angles for generating the inital loop
circle_angle = np.linspace(0, 2*np.pi, 100)
circ1 = r*np.cos(circle_angle)
circ2 = r*np.sin(circle_angle)

e1 = np.zeros((np.size(psi_range), 3))
e2 = np.zeros((np.size(psi_range), 3))
e3 = np.zeros((np.size(psi_range), 3))
xcirc = np.zeros((np.size(psi_range), 100))
ycirc = np.zeros((np.size(psi_range), 100))
zcirc = np.zeros((np.size(psi_range), 100))

for i in range(np.size(psi_range)):
    R1 = np.matrix([[np.cos(psi_range[i]),  np.sin(psi_range[i]), 0],
                   [-np.sin(psi_range[i]), np.cos(psi_range[i]), 0],
                   [0, 0, 1]])
    R2 = np.matrix([[1, 0, 0],
                   [0, np.cos(theta_range[i]), np.sin(theta_range[i])],
                   [0, np.sin(theta_range[i]), np.cos(theta_range[i])]])
    R3 = np.matrix([[np.cos(phi_range[i]),  np.sin(phi_range[i]), 0],
                   [-np.sin(phi_range[i]), np.cos(phi_range[i]), 0],
                   [0, 0, 1]])

    e1[i] = (np.matrix([1, 0, 0])*(R3*R2*R1))
    e2[i] = (np.matrix([0, 1, 0])*(R3*R2*R1))
    e3[i] = (np.matrix([0, 0, 1])*(R3*R2*R1))

    xcirc[i] = x_range[i] + circ1*e1[i, 0] + circ2*e2[i, 0]
    ycirc[i] = y_range[i] + circ1*e1[i, 1] + circ2*e2[i, 1]
    zcirc[i] = z_range[i] + circ1*e1[i ,2] + circ2*e2[i, 2]

disk, = ax.plot(
    xcirc[500], ycirc[500], zcirc[500], color='black', lw=1)  # Initial disk plot

def animate_3D(i):
    disk.set_data(xcirc[i], ycirc[i])  # Update loop data
    disk.set_3d_properties(zcirc[i], 'z')

    return disk,


anim_3D = animation.FuncAnimation(
    fig, animate_3D, frames=int(len(t_range)), interval=t_step, blit=True)

plt.show()
# Uncomment to save animation
'''start = time.time()
anim_3D.save(r'animation/animation_3D.mp4')
end = time.time()
print('3D saving took {} s'.format(round(end-start, 2)))'''
plt.close()