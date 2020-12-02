import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
import time

# Setting parameters
m = 0.5164 # Mass of the pendulum in kilogram
I = 1.45 # Moment of inertia of the pendulum in kilogram per meter^2
k = 2.8 # Spring constant in newton per meter
delta = 7.86 # Torsion constant in newton meter
epsilon = -0.1 # Coupling constant in newton

omega_z = np.sqrt(k/m)
omega_theta = np.sqrt(delta/I)

a = omega_z**2 + omega_theta**2
b = omega_z**2 - omega_theta**2
c = 4*(epsilon**2)/(m*I)

omega_1 = 0.5 * (a + np.sqrt(b**2 + c))
omega_2 = 0.5 * (a - np.sqrt(b**2 + c))

z_0 = 0.01 # Initial vertical displacement in meter
theta_0 = 0.1 # Initial angular displacement in radian

# Generating data
t_range = np.linspace(0, 30, 1000)
B = ((epsilon*z_0)/I + (omega_2**2 - omega_theta**2)*theta_0)/(omega_2**2 - omega_1**2)
D = ((epsilon*z_0)/I + (omega_1**2 - omega_theta**2)*theta_0)/(omega_1**2 - omega_2**2)
theta_range = B*np.cos(omega_1*t_range) + D*np.cos(omega_2*t_range)
z_range = ((delta - I*(omega_1**2))*B*np.cos(omega_1*t_range) + (delta - I*(omega_2**2))*D*np.cos(omega_2*t_range))/epsilon

# Uncomment to use formula given in a paper
'''epsilon = epsilon * (-2)
K = 1 / (omega_1**2-omega_2**2)
theta_range = ((epsilon*z_0)/(2*I))*K*(np.cos(omega_1*t_range) - np.cos(omega_2*t_range)) + theta_0*K*((omega_1**2 - omega_theta**2)*np.cos(omega_2*t_range) - (omega_2**2 - omega_theta**2)*np.cos(omega_1*t_range))
z_range = z_0*k*((omega_1**2 - omega_theta**2)*np.cos(omega_1*t_range) - (omega_2**2 - omega_theta**2)*np.cos(omega_2*t_range)) - (2*I*theta_0/epsilon)*K*(omega_1**2 - omega_theta**2)*(omega_2**2 - omega_theta**2)*(np.cos(omega_1*t_range) - np.cos(omega_2*t_range))'''

# Visualisation separated
fig, axs = plt.subplots(1, 2, figsize=(12, 7))

axs[0].plot(t_range, theta_range)
axs[0].set_xlabel('t/s')
axs[0].set_ylabel('$\theta$/rad')
axs[1].plot(t_range, z_range)
axs[1].set_xlabel('t/s')
axs[1].set_ylabel('$z$/m')

plt.text(
    1.01, 0.05, '$m$ = {} kg'.format(m, 2), transform=plt.gca().transAxes)  # m text
plt.text(
    1.01, 0.00, '$I$ = {} kg/m^2'.format(round(I, 2)), transform=plt.gca().transAxes)  # I text
plt.text(
    1.01, -0.05, '$z_0$ = {} m'.format(round(z_0, 2)), transform=plt.gca().transAxes)  # z_0 text
plt.text(
    1.01, -0.10, '$\\theta_0$ = {} rad'.format(round(theta_0, 2)), transform=plt.gca().transAxes)  # theta_0 text

plt.show()
plt.close()

# Visualisation combined
fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('$t$/s')
ax1.set_ylabel('$\\theta$/rad', color=color)
ax1.plot(t_range, theta_range, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('$z$/m', color=color)
ax2.plot(t_range, z_range, color=color)
ax2.tick_params(axis='y', labelcolor=color)

plt.show()
plt.close()

# 3D animation
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.axes.set_xlim3d(left=-0.15, right=0.15)
ax.axes.set_ylim3d(bottom=-0.15, top=0.15)
ax.axes.set_zlim3d(bottom=1.25*(np.min(z_range)), top=1.25*(np.max(z_range)))
# ax.view_init(azim=55, elev=45) # Set viewing angle

x = np.zeros(10)
y = np.linspace(-0.1, 0.1, 10)
z = np.zeros(10)

line, = ax.plot(x, y, z, color='black', lw=1.5)  # Initial line plot

time_text = ax.text2D(
    0.05, 0.95, "time", transform=ax.transAxes)  # Time text
theta_text = ax.text2D(
    0.05, 0.90, "theta", transform=ax.transAxes)  # theta text
z_text = ax.text2D(
    0.05, 0.85, 'z', transform=ax.transAxes) # z text

def animate_3D(i):
    x_new = y*np.sin(theta_range[i])
    y_new = y*np.cos(theta_range[i])
    z_new = z + z_range[i]
    line.set_data(x_new, y_new)  # Update line data
    line.set_3d_properties(z_new, 'z')

    # Update time text
    time_text.set_text('t = {} s'.format(round(t_range[i], 2)))  
    # Update theta text
    theta_text.set_text('$\\theta$ = {} rad'.format(round(theta_range[i], 2)))
    # Update z text
    z_text.set_text('$z$ = {} m'.format(round(z_range[i], 2)))
    return line, time_text, theta_text, z_text

anim_3D = animation.FuncAnimation(fig, animate_3D, frames=range(int(len(t_range))), interval=10, blit=True)

plt.show()
# Uncomment to save animation
'''start = time.time()
anim_3D.save(r'animation/pendulum.mp4')
end = time.time()
print('3D saving took {} s'.format(round(end-start, 2)))'''
plt.close()