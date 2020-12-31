import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from sympy.solvers import solve
from sympy.solvers.solveset import linsolve
from sympy import *

m = 0.1109862745 # Mass of the pendulum in kilogram
I = 1.058/10000 # Moment of inertia of the pendulum in kilogram per meter^2
k = 2.14 # Spring constant in newton per meter
delta = (2.04)/1000 # Torsion constant in newton meter
epsilon = (9.27/1000)/-2 # Coupling constant in newton
alpha = 0.01 # Friction constant in vertical oscillation
beta = 0.01 # Friction constant in angular oscillation

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
for sol in sols:
    print(sol)

r1 = sols[0]
r2 = sols[1]
r3 = sols[2]
r4 = sols[3]

c = []

for sol in sols:
    ci = m*(sol**2 + 2*lz*sol + omega_z**2)/epsilon
    c.append(ci)

A, B, C, D = symbols('A, B, C, D')

solutions = linsolve(Matrix(([1, 1, 1, 1, theta_0],         [r1, r2, r3, r4, 0], 
                             [c[0], c[1], c[2], c[3], z_0], [c[0]*r1, c[1]*r2, c[2]*r3, c[3]*r4, 0])), 
                             (A, B, C, D))

solutions = list(solutions)[0]

sols = [complex(item) for item in sols]

solutions = [complex(item) for item in solutions]

r1 = sols[0]
r2 = sols[1]
r3 = sols[2]
r4 = sols[3]
