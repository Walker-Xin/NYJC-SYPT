import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

# Set constants
omega = np.pi/0.5
g = 9.81
r = 0.1

A = omega**2
B = g / r

phi_0 = 0
phi_dot_0 = 0.01


def dy_dt(y, t):
    # Solving the equation phi_dotdot = A*sin(phi)*cos(phi) - B*sin(phi)
    # Let y be a vector y = [phi_1, phi_2], where phi_1 = phi and phi_2 = phi_dot
    return [y[1], A*np.sin(y[0])*np.cos(y[0]) - B*np.sin(y[0])]


y_0 = [phi_0, phi_dot_0] # Set initial values
t_range = np.linspace(0, 10, 300)
y_range = integrate.odeint(dy_dt, y_0, t_range)
phi_range = y_range[:,0]

'''print(phi_range)
print(t_range)
print(len(phi_range), len(t_range))'''

plt.xlabel('t/s')
plt.ylabel('$\phi$ /rad')

plt.text(0.85, 0.15, r'$\dot \phi_0$ = {} rad'.format(phi_dot_0), {'fontsize': 10}, transform=plt.gca().transAxes)
plt.text(0.85, 0.10, r'$\omega$ = {} rad/s'.format(round(omega, 2)), {'fontsize': 10}, transform=plt.gca().transAxes)
plt.text(0.85, 0.05, r'$r$ = {} m'.format(round(r, 2)), {'fontsize': 10}, transform=plt.gca().transAxes)

plt.plot(t_range, phi_range)
plt.show()