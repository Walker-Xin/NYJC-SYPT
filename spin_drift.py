import numpy as np
import scipy.integrate as integrate
import scipy.optimize as optimize
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
import time
import sys

# Setting constants
m = 0.1 # Mass of the ring in kilograms
g = 9.81 # Gravitational accelertaion in meter per second^2
r = 0.05 # Radius of the ring in meters
k = 0.2 # Curvature of the parabolic surface

# Initial conditions
D_psi_0 = 0.01
D_theta_0 = 0
D_phi_0 = np.pi
psi_0 = np.radians(0)
theta_0 = np.radians(89)
phi_0 = np.radians(0)
x_0 = 0
y_0 = 0
z_0 = k*(x_0**2 +y_0**2)

'''D_psi_0 = -0.3*np.pi
theta_0 = np.radians(60)
D_phi_0 = -((3/2)*np.sin(theta_0)*(D_psi_0**2) +g/r)/(2*np.tan(theta_0)*(D_psi_0)) # Set conditions to achieve stable state motion'''

t_max=100 # Simulation time in seconds
iterations=10000 # Total number of iterations
t_step=t_max/iterations # Simulation time step
print('Producing simulation with {}s between frames for {} seconds...'.format(t_step, t_max))
t_range = np.linspace(0, t_max, iterations)

# Generating data with numerical DE
# Defining differential equations
def R_solve(R, d):
    return (0.5*k*R*np.sqrt(4*(k**2)*(R**2) + 1) + 0.25*np.log(2*k*R + np.sqrt(4*(k**2)*(R**2) + 1)))/k - d
    

def dy_dt_parabola(w, t):
    D_psi, D_theta, D_phi, psi, theta, phi, x, y, z = w

    DD_psi = (D_theta*D_phi*np.cos(theta) + D_theta*D_phi)/(0.5*np.sin(theta))
    DD_theta = (-(3/2)*(D_psi**2)*np.sin(theta)*np.cos(theta) - 2*D_phi*D_psi*np.sin(theta) - (g/r)*np.cos(theta))/(1.5)
    DD_phi = (3*D_theta*D_psi*np.sin(theta) - 2*np.cos(theta)*DD_psi)/(2)
    D_x_ = -r*(D_psi*np.cos(psi)*np.cos(theta) - D_theta*np.sin(psi)*np.sin(theta) + D_phi*np.cos(psi))
    D_y_ = -r*(D_psi*np.sin(psi)*np.cos(theta) + D_theta*np.cos(psi)*np.sin(theta) + D_phi*np.sin(psi))

    alpha = np.arctan2(y, x)
    R = np.sqrt(z/k)

    D_R = D_y_*np.sin(alpha) + D_x_*np.cos(alpha)
    D_t = D_y_*np.cos(alpha) - D_x_*np.sin(alpha)

    D_x = -D_t*np.sin(alpha)
    D_y = D_t*np.cos(alpha)
    D_z = 2*k*R*D_R

    dy_dt = [DD_psi, DD_theta, DD_phi, D_psi, D_theta, D_phi, D_x, D_y, D_z]

    return dy_dt


def dy_dt(w, t):
    D_psi, D_theta, D_phi, psi, theta, phi, x, y, z = w

    DD_psi = (D_theta*D_phi*np.cos(theta) + D_theta*D_phi)/(0.5*np.sin(theta))
    DD_theta = (-(3/2)*(D_psi**2)*np.sin(theta)*np.cos(theta) - 2*D_phi*D_psi*np.sin(theta) - (g/r)*np.cos(theta))/(1.5)
    DD_phi = (3*D_theta*D_psi*np.sin(theta) - 2*np.cos(theta)*DD_psi)/(2)
    D_x = -r*(D_psi*np.cos(psi)*np.cos(theta) - D_theta*np.sin(psi)*np.sin(theta) + D_phi*np.cos(psi))
    D_y = -r*(D_psi*np.sin(psi)*np.cos(theta) + D_theta*np.cos(psi)*np.sin(theta) + D_phi*np.sin(psi))

    dy_dt = [DD_psi, DD_theta, DD_phi, D_psi, D_theta, D_phi, D_x, D_y, 0]

    return dy_dt


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

start = time.time()
w_0 = [D_psi_0, D_theta_0, D_phi_0, psi_0, theta_0, phi_0, x_0, y_0, z_0] # Set initial conditions

w_range = integrate.odeint(dy_dt, w_0, t_range) # Solving equation with ODEint solver

psi_dot_range = w_range[:, 0]
theta_dot_range = w_range[:, 1]
phi_dot_range = w_range[:, 2]
psi_range = w_range[:, 3]
theta_range = w_range[:, 4]
phi_range = w_range[:, 5]
x_range = w_range[:, 6]
y_range = w_range[:, 7]
z_range = w_range[:, 8]

end = time.time()
print('Data generation took {} s'.format(round(end-start, 2)))

# Motion visualisation
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
ax.axes.set_xlim3d(left=np.min(x_range)-r, right=np.max(x_range)+r)
ax.axes.set_ylim3d(bottom=np.min(y_range)-r, top=np.max(y_range)+r)
ax.axes.set_zlim3d(bottom=-5*(2*r), top=5*(2*r))
ax.set_box_aspect((1, 1, 1), zoom=1)
# ax.view_init(azim=0, elev=45) # Set viewing angle

# A set of angles for generating the inital loop
circle_angle = np.linspace(0, 2*np.pi, 100)
circ1 = r*np.cos(circle_angle)
circ2 = r*np.sin(circle_angle)

# Initialise data sets
e2pp = np.zeros((np.size(psi_range), 3))
e1 = np.zeros((np.size(psi_range), 3))
e2 = np.zeros((np.size(psi_range), 3))
e3 = np.zeros((np.size(psi_range), 3))
xcirc = np.zeros((np.size(psi_range), 100))
ycirc = np.zeros((np.size(psi_range), 100))
zcirc = np.zeros((np.size(psi_range), 100))
x_path = np.zeros(np.size(psi_range))
y_path = np.zeros(np.size(psi_range))
z_path = np.zeros(np.size(psi_range))
x_A = np.zeros(np.size(psi_range))
y_A = np.zeros(np.size(psi_range))
z_A = np.zeros(np.size(psi_range))

for i in range(np.size(psi_range)):
    # Setting Euler angles
    R1 = np.matrix([[np.cos(psi_range[i]),  np.sin(psi_range[i]), 0],
                   [-np.sin(psi_range[i]), np.cos(psi_range[i]), 0],
                   [0, 0, 1]])
    R2 = np.matrix([[1, 0, 0],
                   [0, np.cos(theta_range[i]), np.sin(theta_range[i])],
                   [0, np.sin(theta_range[i]), np.cos(theta_range[i])]])
    R3 = np.matrix([[np.cos(phi_range[i]),  np.sin(phi_range[i]), 0],
                   [-np.sin(phi_range[i]), np.cos(phi_range[i]), 0],
                   [0, 0, 1]])

    e2pp[i] = (np.matrix([0, 1, 0])*(R2*R1))

    e1[i] = (np.matrix([1, 0, 0])*(R3*R2*R1))
    e2[i] = (np.matrix([0, 1, 0])*(R3*R2*R1))
    e3[i] = (np.matrix([0, 0, 1])*(R3*R2*R1))

    xcirc[i] = x_range[i] + circ1*e1[i, 0] + circ2*e2[i, 0]
    ycirc[i] = y_range[i] + circ1*e1[i, 1] + circ2*e2[i, 1]
    zcirc[i] = z_range[i] + circ1*e1[i ,2] + circ2*e2[i, 2]

    x_path[i] = x_range[i] - r*e2pp[i, 0]
    y_path[i] = y_range[i] - r*e2pp[i, 1]
    z_path[i] = z_range[i]

    x_A[i] = x_range[i] + r*e2pp[i, 0]
    y_A[i] = y_range[i] + r*e2pp[i, 1]
    z_A[i] = z_range[i] + r*e2pp[i, 2]

disk, = ax.plot(
    xcirc[0], ycirc[0], zcirc[0], color='black', lw=1)  # Initial disk plot
path, = ax.plot(
    x_path[0], y_path[0], z_path[0], color='blue', lw=0.7) # Initial path plot
point, = ax.plot(
    x_A[0], y_A[0], z_A[0], marker='o', color='red') # Initial fixed point plot

def animate_3D(i):
    disk.set_data(xcirc[i], ycirc[i])  # Update loop data
    disk.set_3d_properties(zcirc[i], 'z')

    path.set_data(x_path[:i], y_path[:i]) # Updated path data
    path.set_3d_properties(z_path[:i], 'z')

    point.set_data(x_A[i], y_A[i]) # Updated fixed point data
    point.set_3d_properties(z_A[i], 'z')

    return disk, path, point,


anim_3D = animation.FuncAnimation(
    fig, animate_3D, frames=int(len(t_range)), interval=t_step, blit=True)

plt.show()
# Uncomment to save animation
'''start = time.time()
anim_3D.save(r'animation/animation_3D.mp4')
end = time.time()
print('3D saving took {} s'.format(round(end-start, 2)))'''
plt.close()

# Transformation to parabolic surface
x_p_range = np.zeros(np.size(x_range))
y_p_range = np.zeros(np.size(x_range))
z_p_range = np.zeros(np.size(x_range))

def R_solve(R, d):
    return (0.5*k*R*np.sqrt(4*(k**2)*(R**2) + 1) + 0.25*np.log(2*k*R + np.sqrt(4*(k**2)*(R**2) + 1)))/k - d

for i in range(np.size(x_range)):
    d = np.sqrt(x_range[i]**2 + y_range[i]**2)
    sol = optimize.fsolve(R_solve, 0, args = (d))
    R = sol[0]
    alpha = np.arctan2(y_range[i], x_range[i])

    x_p = R*np.cos(alpha)
    y_p = R*np.sin(alpha)
    z_p = k*(R**2)

    x_p_range[i] = x_p
    y_p_range[i] = y_p
    z_p_range[i] = z_p

fig, axs = plt.subplots(1, 1, figsize=(12, 12))

color = 'tab:red'
axs.plot(x_p_range, z_p_range, color=color)

plt.show()
plt.close()

# 3D animation
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.axes.set_xlim3d(left=-0.5, right=0.5)
ax.axes.set_ylim3d(bottom=-0.5, top=0.5)
ax.axes.set_zlim3d(bottom=0, top=1)
ax.set_box_aspect((1, 1, 1), zoom=1)
# ax.view_init(azim=0, elev=45) # Set viewing angle

# Plot parabolic surface
u = np.linspace(0, 0.5, 100)
v = np.linspace(0, 2*np.pi, 100)
U, V = np.meshgrid(u, v)

Z = k*(U**2)
X = U*np.cos(V)
Y = U*np.sin(V)
ax.plot_surface(
    X, Y, Z, color='gray', alpha=0.8)

# A set of angles for generating the inital loop
circle_angle = np.linspace(0, 2*np.pi, 100)
circ1 = r*np.cos(circle_angle)
circ2 = r*np.sin(circle_angle)

# Initialise data sets
e2pp = np.zeros((np.size(psi_range), 3))
e1 = np.zeros((np.size(psi_range), 3))
e2 = np.zeros((np.size(psi_range), 3))
e3 = np.zeros((np.size(psi_range), 3))
xcirc = np.zeros((np.size(psi_range), 100))
ycirc = np.zeros((np.size(psi_range), 100))
zcirc = np.zeros((np.size(psi_range), 100))
x_path = np.zeros(np.size(psi_range))
y_path = np.zeros(np.size(psi_range))
z_path = np.zeros(np.size(psi_range))
x_A = np.zeros(np.size(psi_range))
y_A = np.zeros(np.size(psi_range))
z_A = np.zeros(np.size(psi_range))

for i in range(np.size(psi_range)):
    # Setting Euler angles
    R1 = np.matrix([[np.cos(psi_range[i]),  np.sin(psi_range[i]), 0],
                   [-np.sin(psi_range[i]), np.cos(psi_range[i]), 0],
                   [0, 0, 1]])
    R2 = np.matrix([[1, 0, 0],
                   [0, np.cos(theta_range[i]), np.sin(theta_range[i])],
                   [0, np.sin(theta_range[i]), np.cos(theta_range[i])]])
    R3 = np.matrix([[np.cos(phi_range[i]),  np.sin(phi_range[i]), 0],
                   [-np.sin(phi_range[i]), np.cos(phi_range[i]), 0],
                   [0, 0, 1]])

    e2pp[i] = (np.matrix([0, 1, 0])*(R2*R1))

    e1[i] = (np.matrix([1, 0, 0])*(R3*R2*R1))
    e2[i] = (np.matrix([0, 1, 0])*(R3*R2*R1))
    e3[i] = (np.matrix([0, 0, 1])*(R3*R2*R1))

    xcirc[i] = x_p_range[i] + circ1*e1[i, 0] + circ2*e2[i, 0]
    ycirc[i] = y_p_range[i] + circ1*e1[i, 1] + circ2*e2[i, 1]
    zcirc[i] = z_p_range[i] + circ1*e1[i ,2] + circ2*e2[i, 2]

    x_path[i] = x_p_range[i] - r*e2pp[i, 0]
    y_path[i] = y_p_range[i] - r*e2pp[i, 1]
    z_path[i] = z_p_range[i]

    x_A[i] = x_p_range[i] + r*e2pp[i, 0]
    y_A[i] = y_p_range[i] + r*e2pp[i, 1]
    z_A[i] = z_p_range[i] + r*e2pp[i, 2]

disk, = ax.plot(
    xcirc[0], ycirc[0], zcirc[0], color='black', lw=1)  # Initial disk plot
'''path, = ax.plot(
    x_path[0], y_path[0], z_path[0], color='blue', lw=0.7) # Initial path plot
point, = ax.plot(
    x_A[0], y_A[0], z_A[0], marker='o', color='red') # Initial fixed point plot'''

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