import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from sympy.solvers import solve
from sympy.solvers.solveset import linsolve
from sympy import *
import time

m = 0.5164 # Mass of the pendulum in kilogram
I = 1.45/10000 # Moment of inertia of the pendulum in kilogram per meter^2
k = 2.8 # Spring constant in newton per meter
delta = (7.86)/10000 # Torsion constant in newton meter
epsilon = (9.27/1000)/-2 # Coupling constant in newton
alpha = 0.01 # Friction constant in vertical oscillation
beta= 0.0001 # Friction constant in angular oscillation

lz = alpha/(2*m)
lt = beta/(2*I)

omega_z = np.sqrt(k/m)
omega_theta = np.sqrt(delta/I)

a = omega_z**2 + omega_theta**2
b = omega_z**2 - omega_theta**2
c = 4*(epsilon**2)/(m*I)

omega_1 = np.sqrt(0.5 * (a + np.sqrt(b**2 + c)))
omega_2 = np.sqrt(0.5 * (a - np.sqrt(b**2 + c)))

z_0 = 0.1 # Initial vertical displacement in meter
theta_0 = 0/180 * np.pi # Initial angular displacement in radian

t_max=100 # simulation time in seconds
iterations=100000 # total number of iterations
t_step=t_max/iterations # simulation time step
print('Producing simulation with {}s between frames...'.format(t_step))
t_range = np.linspace(0, t_max, iterations)

def f(r):
    return r**4 + 2*(lz + lt)*(r**3) + (omega_z**2 + omega_theta**2 + 4*lz*lt)*(r**2) + 2*(lz*omega_theta**2 + lt*omega_z**2)*r + (omega_z*omega_theta)**2 - (epsilon**2)/(m*I)


def f_p(r):
    return 4*(r**3) + 6*(lz + lt)*(r**2) + 2*(omega_z**2 + omega_theta**2 + 4*lz*lt)*(r) + 2*(lz*omega_theta**2 + lt*omega_z**2)


def f_pp(r):
    return 12*(r**2) + 12*(lz + lt)*(r) + 2*(omega_z**2 + omega_theta**2 + 4*lz*lt)


phi_sample = np.linspace(-500, 500, 10000)

y = np.zeros(np.size(phi_sample))

for i in range(len(phi_sample)):
    y[i] = f(phi_sample[i])

'''plt.plot(phi_sample, y)
plt.plot(phi_sample, np.zeros(np.size(phi_sample)))
plt.show()'''

phi = Symbol('phi')
sols = solve(f(phi), phi)

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

print('Eigenfrequencies are: \n {} \n {} \n {} \n {}'.format(r1, r2, r3, r4))

solutions = [complex(item) for item in solutions]

A = solutions[0]
B = solutions[1]
C = solutions[2]
D = solutions[3]

c1 = complex(c1)
c2 = complex(c2)
c3 = complex(c3)
c4 = complex(c4)

t = Symbol('t')

z_range = np.real(A*np.exp(r1*t_range) + B*np.exp(r2*t_range) + C*np.exp(r3*t_range) + D*np.exp(r4*t_range))
theta_range = np.real(c1*A*np.exp(r1*t_range) + c2*B*np.exp(r2*t_range) + c3*C*np.exp(r3*t_range) + c4*D*np.exp(r4*t_range))
z_dot_range = np.real(A*r1*np.exp(r1*t_range) + B*r2*np.exp(r2*t_range) + C*r3*np.exp(r3*t_range) + D*r4*np.exp(r4*t_range))
theta_dot_range = np.real(c1*A*r1*np.exp(r1*t_range) + c2*B*r2*np.exp(r2*t_range) + c3*C*r3*np.exp(r3*t_range) + c4*D*r4*np.exp(r4*t_range))

E_0 = 0.5*(k*(z_0**2) + delta*(theta_0**2))

start = time.time()
import scipy.integrate as integrate

energy_loss_range = alpha*(z_dot_range**2) + beta*(theta_dot_range**2)

for i in range(iterations):
    loss = integrate.trapz(energy_loss_range[:i], dx=t_step)
    if abs((loss-E_0)/E_0) < 0.001:
        t_e = (i/iterations)*t_max
        print('At t = {}s, the oscillations virtually stoped.'.format(round(t_e, 2)))
        break
    else:
        pass

print('Time used: {}s'.format(round(time.time() - start, 2)))

'''trans = quad(transfer, 0, 14)[0]
print(trans)

print('Percentaged transferred = {}%'.format(round((trans/E_0)*100), 2))'''