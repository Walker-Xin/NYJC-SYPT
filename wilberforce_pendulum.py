import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
from sympy.solvers import solve
from sympy.solvers.solveset import linsolve
from sympy import symbols, Symbol, Matrix
import time

# Setting parameters
m = 0.1109862745 # Mass of the pendulum in kilogram
I = 1.058/10000 # Moment of inertia of the pendulum in kilogram per meter^2
k = 2.14 # Spring constant in newton per meter
delta = (2.04)/1000 # Torsion constant in newton meter
epsilon = (9.27/1000)/-2 # Coupling constant in newton
alpha = 0.008 # Friction constant in vertical oscillation
beta= 0.00008 # Friction constant in angular oscillation

lz = alpha/(2*m)
lt = beta/(2*I)

omega_z = np.sqrt(k/m)
omega_theta = np.sqrt(delta/I)

a = omega_z**2 + omega_theta**2
b = omega_z**2 - omega_theta**2
c = 4*(epsilon**2)/(m*I)

omega_1 = np.sqrt(0.5 * (a + np.sqrt(b**2 + c)))
omega_2 = np.sqrt(0.5 * (a - np.sqrt(b**2 + c)))

z_0 = 0.01 # Initial vertical displacement in meter
theta_0 = 30/180 * np.pi # Initial angular displacement in radian

t_max=100 # simulation time in seconds
iterations=10000 # total number of iterations
t_step=t_max/iterations # simulation time step
print('Producing simulation with {}s between frames...'.format(t_step))
t_range = np.linspace(0, t_max, iterations)

# Generating data with Euler's method
'''theta_range=np.zeros(iterations)
theta_range[0]=theta_0
z_range=np.zeros(iterations)
z_range[0]=z_0

theta_dot_range=np.zeros(iterations)
theta_dot_range[0]=0
z_dot_range=np.zeros(iterations)
z_dot_range[0]=0

theta_dotdot_range=np.zeros(iterations)
theta_dotdot_range[0]=(epsilon*z_0-delta*theta_0)/I
z_dotdot_range=np.zeros(iterations)
z_dotdot_range[0]=(epsilon*theta_0-k*z_0)/m

for n in range(1,iterations):
    theta_range[n]=theta_range[n-1]+theta_dot_range[n-1]*t_step # theta
    z_range[n]=z_range[n-1]+z_dot_range[n-1]*t_step # z

    theta_dot_range[n]=theta_dot_range[n-1]+theta_dotdot_range[n-1]*t_step # theta_dot
    z_dot_range[n]=z_dot_range[n-1]+z_dotdot_range[n-1]*t_step # z_dot

    theta_dotdot_range[n]=(epsilon*z_range[n]-delta*theta_range[n])/I # theta_dotdot
    z_dotdot_range[n]=(epsilon*theta_range[n]-k*z_range[n])/m #z_dotdot'''

# Generating data with numerical DE
'''def dy_dt(y, t): # Defining differential equations
    z_1, z_2, theta_1, theta_2 = y

    dy_dt = [z_2, (epsilon/m)*theta_1-(k/m)*z_1-alpha*z_2, theta_2, (epsilon/I)*z_1-(delta/I)*theta_1-beta*theta_2]
    return dy_dt


y_0 = [z_0, 0, theta_0, 0] # Setting intial conditions
y_range = integrate.odeint(dy_dt, y_0, t_range)
z_range = y_range[:, 0]
z_dot_range = y_range[:, 1]
theta_range = y_range[:, 2]
theta_dot_range = y_range[:, 3]'''

# Generating data with analytical solution
'''B = ((epsilon*z_0)/I + (omega_2**2 - omega_theta**2)*theta_0)/(omega_2**2 - omega_1**2)
D = ((epsilon*z_0)/I + (omega_1**2 - omega_theta**2)*theta_0)/(omega_1**2 - omega_2**2)
theta_range = B*np.cos(omega_1*t_range) + D*np.cos(omega_2*t_range)
z_range = ((delta - I*(omega_1**2))*B*np.cos(omega_1*t_range) + (delta - I*(omega_2**2))*D*np.cos(omega_2*t_range))/epsilon
theta_dot_range = -B*omega_1*np.cos(omega_1*t_range) - D*omega_2*np.cos(omega_2*t_range)
z_dot_range = (-(delta - I*(omega_1**2))*B*omega_1*np.cos(omega_1*t_range) - (delta - I*(omega_2**2))*D*omega_2*np.cos(omega_2*t_range))/epsilon'''

# Generating data with damping with analytical solution
def f(r):
    return r**4 + 2*(lz + lt)*(r**3) + (omega_z**2 + omega_theta**2 + 4*lz*lt)*(r**2) + 2*(lz*omega_theta**2 + lt*omega_z**2)*r + (omega_z*omega_theta)**2 - (epsilon**2)/(m*I)

phi = Symbol('phi')
sols = solve(f(phi), phi)
for sol in sols:
    print(sol)

r1 = sols[0]
r2 = sols[1]
r3 = sols[2]
r4 = sols[3]

c1 = m*(r1**2 + 2*lz*r1 + omega_z**2)/epsilon
c2 = m*(r2**2 + 2*lz*r2 + omega_z**2)/epsilon
c3 = m*(r3**2 + 2*lz*r3 + omega_z**2)/epsilon
c4 = m*(r4**2 + 2*lz*r4 + omega_z**2)/epsilon

A, B, C, D = symbols('A, B, C, D')

solutions = linsolve(Matrix(([1, 1, 1, 1, z_0],         [r1, r2, r3, r4, 0], 
                             [c1, c2, c3, c4, theta_0], [c1*r1, c2*r2, c3*r3, c4*r4, 0])), 
                             (A, B, C, D))

solutions = list(solutions)[0]

# Convert Sympy complex format to Python built-in complex format
r1 = complex(r1)
r2 = complex(r2)
r3 = complex(r3)
r4 = complex(r4)

solutions = [complex(item) for item in solutions]

A = solutions[0]
B = solutions[1]
C = solutions[2]
D = solutions[3]

c1 = complex(c1)
c2 = complex(c2)
c3 = complex(c3)
c4 = complex(c4)

z_range = np.real(A*np.exp(r1*t_range) + B*np.exp(r2*t_range) + C*np.exp(r3*t_range) + D*np.exp(r4*t_range))
theta_range = np.real(c1*A*np.exp(r1*t_range) + c2*B*np.exp(r2*t_range) + c3*C*np.exp(r3*t_range) + c4*D*np.exp(r4*t_range))
z_dot_range = np.real(A*r1*np.exp(r1*t_range) + B*r2*np.exp(r2*t_range) + C*r3*np.exp(r3*t_range) + D*r4*np.exp(r4*t_range))
theta_dot_range = np.real(c1*A*r1*np.exp(r1*t_range) + c2*B*r2*np.exp(r2*t_range) + c3*C*r3*np.exp(r3*t_range) + c4*D*r4*np.exp(r4*t_range))

# Visualisation separated
fig, axs = plt.subplots(2, 2, figsize=(12, 7))

color = 'tab:red'
axs[0][0].plot(t_range, theta_range, color=color)
axs[0][0].set_xlabel('t/s')
axs[0][0].set_ylabel('$\\theta$/rad')
color = 'coral'
axs[1][0].plot(t_range, theta_dot_range, color=color)
axs[1][0].set_xlabel('t/s')
axs[1][0].set_ylabel('$\dot \\theta$/rad$\cdot$s$^-1$')

color = 'tab:blue'
axs[0][1].plot(t_range, 100*z_range, color=color)
axs[0][1].set_xlabel('t/s')
axs[0][1].set_ylabel('$z$/cm')
color = 'lightskyblue'
axs[1][1].plot(t_range, 100*z_dot_range, color=color)
axs[1][1].set_xlabel('t/s')
axs[1][1].set_ylabel('$\dot z$/cm$\cdot$s$^-1$')

plt.text(
    1.01, 0.075, '$m$ = {} kg'.format(m, 2), transform=plt.gca().transAxes)  # m text
plt.text(
    1.01, 0.00, '$I$ = {} kg/m$^2$'.format(round(I, 2)), transform=plt.gca().transAxes)  # I text
plt.text(
    1.01, -0.075, '$z_0$ = {} cm'.format(round(100*z_0, 2)), transform=plt.gca().transAxes)  # z_0 text
plt.text(
    1.01, -0.15, '$\\theta_0$ = {} rad'.format(round(theta_0, 2)), transform=plt.gca().transAxes)  # theta_0 text

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
ax2.set_ylabel('$z$/cm', color=color)
ax2.plot(t_range, 100*z_range, color=color)
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
# ax.view_init(azim=0, elev=90) # Set viewing angle

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
    theta_text.set_text('$\\theta$ = {} rad = {} degrees'.format(round(theta_range[i], 2), round(np.degrees(theta_range[i]), 1)))
    # Update z text
    z_text.set_text('$z$ = {} cm'.format(round(100*z_range[i], 2)))
    return line, time_text, theta_text, z_text

anim_3D = animation.FuncAnimation(fig, animate_3D, frames=range(int(len(t_range))), interval=t_step*1000, blit=True)

plt.show()
# Uncomment to save animation
start = time.time()
'''anim_3D.save(r'animation/pendulum.mp4')
end = time.time()
print('3D saving took {} s'.format(round(end-start, 2)))'''
plt.close()