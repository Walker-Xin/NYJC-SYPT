import numpy as np
import scipy.integrate as integrate
import scipy.optimize as optimize
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
import time
import csv
import os

# Setting constants
RPM = 151
omega = RPM/60*2*np.pi # Angular velocity of the loop in radian per second
g = 9.81  # Gravitational accelertaion in meter per second^2
r = 7.675/100  # Radius of the loop in meter
R = 0.9/100 # Radius of the bead in meter
k = 1 + 2*(r**2)/(5*((r-R)**2)) # Correction factor
m = 24.17/1000 # Mass of the bead in kilogram
b = 0.065 # Frcition constant 

A = omega**2
if R == 0:
    B = g/r
else:
    B = g/(r-R)

phi_0 = 0  # Initial angular displacement from -y-axis anticlockwise of the bead in radian
phi_dot_0 = 0.01  # Initial angular velocity from -y-axis anticlockwise of the bead in radian per second

t_max=10 # Simulation time in seconds
iterations=15000 # Total number of iterations
t_step=t_max/iterations # Simulation time step
print('Producing simulation with {}s between frames for {} seconds...'.format(t_step, t_max))
t_range = np.linspace(0, t_max, iterations)

# Generating data with numerical DE
# Defining differential equations
def wet_friction_ball(y, t):
    # Solving the equation of motion: phi_dotdot = A*sin(phi)*cos(phi) - B*sin(phi) - k*phi_dot
    return [y[1], (A*np.sin(y[0])*np.cos(y[0]) - B*np.sin(y[0]) - (b/m)*y[1])/k]


start = time.time()
y_0 = [phi_0, phi_dot_0]  # Setting initial conditions

y_range = integrate.odeint(wet_friction_ball, y_0, t_range) # Solving equation with ODEint solver
phi_range = y_range[:, 0]
phi_dot_range = y_range[:, 1]

theta_range = -omega * t_range

# Compute effective energy
def E_equiv(p, pd):
    return 0.5*m*((r-R)**2)*(pd**2)*k - 0.5*m*((r-R)**2)*(omega**2)*(np.sin(p)**2) + m*g*(r-R)*(1-np.cos(p))


def U_equiv(p, pd):
    return - 0.5*m*((r-R)**2)*(omega**2)*(np.sin(p)**2) + m*g*(r-R)*(1-np.cos(p))


E_range = E_equiv(phi_range, phi_dot_range)
U_range = U_equiv(phi_range, phi_dot_range)

# Calculate projected rest position and time
if B/A < 1:
    phi_e = np.arccos(B/A)
    print('Equilibrium angle is {}rad'.format(round(phi_e, 3)))
    energy_loss_range = -b*np.square((r-R)*phi_dot_range)

    E_0 = 0.5*m*((r-R)**2)*(phi_dot_0**2)*k - 0.5*m*((r-R)**2)*(omega**2)*(np.sin(phi_0)**2) + m*g*(r-R)*(1-np.cos(phi_0))

    E_f = -0.5*m*((r-R)**2)*(omega**2)*(np.sin(phi_e)**2) + m*g*(r-R)*(1-np.cos(phi_e))

    delta = (E_0 - E_f) # Loss of energy

    for i in range(iterations):
        loss = integrate.trapz(energy_loss_range[:i], dx=t_step)
        ratio = 1 - abs(loss/delta)
        if ratio < 0.0001:
            t_e = (i/iterations)*t_max
            break
    else:
        pass
else:
    print('The bead is unable to swing up')
    phi_e = 0

if 't_e' in globals():
    print('Around t = {}s, the oscillations virtually stoped.'.format(round(t_e, 2)))
else:
    print('t_e not found!')
    t_e = 0

end = time.time()
print('Data generation took {} s'.format(round(end-start, 2)))

# Motion visualisation
fig, axs = plt.subplots(1, 2, figsize=(12, 7))

axs[0].plot(t_range, phi_range)
axs[0].set_xlabel('t/s')
axs[0].set_ylabel('$\phi$/rad')
axs[1].plot(t_range, phi_dot_range)
axs[1].set_xlabel('t/s')
axs[1].set_ylabel('$\dot \phi$/rad$\cdot$s$^-1$')

# Add projected rest position and time
phi_e_range = np.full(np.shape(t_range), phi_e) 
color = 'tab:red'
axs[0].plot(t_range, phi_e_range, color=color, linestyle=':') # Projected rest position
axs[0].plot(t_e, phi_e, color=color, marker='x')
axs[1].plot(t_e, 0, color=color, marker='x')

plt.text(
    1.01, 0.05, r'$\dot \phi_0$ = {} rad$\cdot$s$^-1$'.format(phi_dot_0), transform=plt.gca().transAxes)  # phi_dot_0 text
plt.text(
    1.01, 0.00, r'$\omega$ = {} rad$\cdot$s$^-1$'.format(round(omega, 2)), transform=plt.gca().transAxes)  # theta_dot_0 text
plt.text(
    1.01, -0.05, r'$r$ = {} m'.format(round(r, 2)), transform=plt.gca().transAxes)  # Radius text

plt.show()
plt.close()

# Phase diagram
omega_range = np.linspace(0.01, np.pi/0.15, 1000)
phi_e_range_plus = np.zeros(np.shape(omega_range))
phi_e_range_minus = np.zeros(np.shape(omega_range))
for i in range(np.size(omega_range)):
    a = omega_range[i]**2
    if a > B:
        phi_e_range_plus[i] = np.arccos(B/a)
        phi_e_range_minus[i] = -np.arccos(B/a)
    else:
        pass

fig, axs = plt.subplots(1, 1, figsize=(12, 7))

axs.plot(omega_range, phi_e_range_minus, color='tab:blue')
axs.plot(omega_range, phi_e_range_plus, color='tab:blue')
axs.plot(np.sqrt(B), 0, color='tab:red', marker='x')
axs.set_xlabel('$\omega$/rad$\cdot$s$^-1$')
axs.set_ylabel('$\phi_e$/rad')

plt.show()
plt.close()

# Energy Visualisation
fig, axs = plt.subplots(1, 1, figsize=(8, 6))

axs.plot(t_range, E_range, color = 'tab:blue', label='E')
axs.plot(t_range, U_range, color = 'tab:red', label='U')
axs.plot(t_range, np.full(np.shape(E_range), E_f), color = 'tab:gray')
index = int(t_e/t_range.max()*len(t_range))
axs.plot(t_e, E_range[index], color = 'black', marker='x')
axs.set_xlabel('t/s')
axs.set_ylabel('$E_{eff}$/J')
axs.legend()
plt.text(0.85, 0.01, r'$r$ = {} cm'.format(round(r*100, 2)), transform=plt.gca().transAxes)  # Radius text
plt.text(0.85, 0.06, r'$R$ = {} cm'.format(round(R*100, 2)), transform=plt.gca().transAxes)  # Radius text

plt.show()
plt.close()

import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3

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
phi_text = ax.text(
    0.05, 0.90, "phi", transform=ax.transAxes)  # phi text
phi_dot_text = ax.text(
    0.05, 0.85, "phi_dot", transform=ax.transAxes)  # phi_dot text
theta_dot_text = ax.text(
    0.05, 0.80, "$\dot \\theta$ = $\omega$ = {} rad$\cdot$s$^-1$".format(round(omega, 2)), transform=ax.transAxes)  # theta_dot text


def animate_2D(i):
    x_coord = x[i]
    y_coord = y[i]
    point.set_data([x_coord], [y_coord])  # Update bead data

    # Update time text
    time_text.set_text(
        't = {} s'.format(round(t_range[i], 2)))
    # Update phi text
    phi_text.set_text(
        '$\phi$ = {} rad = {} degrees'.format(round(phi_range[i], 2), round(np.degrees(phi_range[i]), 1)))
    # Update phi_dot text
    phi_dot_text.set_text(
        '$\dot \phi$ = {} rad$\cdot$s$^-1$'.format(round(phi_dot_range[i], 2)))
    return point, time_text, phi_text, phi_dot_text


anim_2D = animation.FuncAnimation(
    fig, animate_2D, frames=int(len(t_range)), interval=t_step, blit=True)

plt.show()
# Uncomment to save animation
'''start = time.time()
anim_2D.save(r'animation/animation_2D.mp4')
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
ax.view_init(azim=0, elev=90) # Set viewing angle

# A set of angles for generating the inital loop
circle_angle = np.linspace(0, 2*np.pi, 100)

x = r * np.sin(phi_range) * np.sin(theta_range)
y = r * np.sin(phi_range)
z = -r * np.cos(phi_range)

point, = ax.plot([0], [0], [0], marker='.', color='r')  # Initial bead plot

x_circle = np.zeros(np.shape(circle_angle))
y_circle = r * np.cos(circle_angle)
z_circle = r * np.sin(circle_angle)

circle, = ax.plot(
    x_circle, y_circle, z_circle, color='black', lw=0.5)  # Initial loop plot

time_text = ax.text2D(
    0.05, 0.95, "Time", transform=ax.transAxes)  # Time text
phi_text = ax.text2D(
    0.05, 0.90, "phi", transform=ax.transAxes)  # phi text
phi_dot_text = ax.text2D(
    0.05, 0.85, "phi_dot", transform=ax.transAxes)  # phi_dot text
theta_dot_text = ax.text2D(
    0.05, 0.80, "$\dot \\theta$ = $\omega$ = {} rad$\cdot$s$^-1$".format(round(omega, 2)), transform=ax.transAxes)  # theta_dot text


def animate_3D(i):
    x_coord = x[i]
    y_coord = y[i]
    z_coord = z[i]
    point.set_data([0], [y_coord])  # Update bead data
    point.set_3d_properties([z_coord], 'z')

    x_circle_new = y_circle * np.sin(theta_range[i])
    y_circle_new = y_circle * np.cos(theta_range[i])
    z_circle_new = z_circle
    circle.set_data(x_circle_new, y_circle_new)  # Update loop data
    circle.set_3d_properties(z_circle_new, 'z')

    # Update time text
    time_text.set_text('t = {} s'.format(
        round(t_range[i], 2)))
    # Update phi text
    phi_text.set_text(
        '$\phi$ = {} rad = {} degrees'.format(round(phi_range[i], 2), round(np.degrees(phi_range[i]), 1)))
    # Update phi_dot text
    phi_dot_text.set_text(
        '$\dot \phi$ = {} rad$\cdot$s$^-1$'.format(round(phi_dot_range[i], 2)))
    return point, circle, time_text, phi_text, phi_dot_text, ax


anim_3D = animation.FuncAnimation(
    fig, animate_3D, frames=int(len(t_range)), interval=t_step, blit=True)

plt.show()
# Uncomment to save animation
'''start = time.time()
anim_3D.save(r'animation/animation_3D.mp4')
end = time.time()
print('3D saving took {} s'.format(round(end-start, 2)))'''
plt.close()