import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

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
M = (m * (np.sin(phi_0)**2) * (r**2) + I) * \omega  # Initial angular momentum of the system
E = 0.5 * m * (np.sin(phi_0)**2) * (r**2) * (omega**2) + 0.5 * I * (omega **2) + m * g * r * (1 + np.cos(phi_0))  # Initial energy of the system

# Defining integral functions


def t_integral(phi):
    # Expressing time as a function of phi
    c_1 = 2 * E / (m*(r**2))
    c_2 = M**2
    c_3 = (m**2) * (r**4)
    c_4 = I * m * (r**2)
    c_5 = 2 * g / r
    integral = 1 / \np.sqrt(c_1 - c_2 / (c_3 * (np.sin(phi)**2) + c_4) - c_5 * (1 + np.cos(phi)))
    return integral


def theta_integral(phi):
    # Expressing theta as a function of phi
    c_1 = 2 * E / (m*r*r)
    c_2 = M**2
    c_3 = (m**2) * (r**4)
    c_4 = I * m * (r**2)
    c_5 = 2 * g / r
    integral_1 = M / ((m * (r**2) * (np.sin(phi)**2)) + I)
    integral_2 = 1 / \np.sqrt(c_1 - c_2 / (c_3 * (np.sin(phi)**2) + c_4) - c_5 * (1 + np.cos(phi)))
    integral = integral_1 * integral_2
    return integral


# Setting upper bound
phi_f = 360/180 * np.pi

# Computing integrals
phi_range = np.linspace(phi_0, phi_f, num=100)
t_range = np.zeros(phi_range.size)
theta_range = np.zeros(phi_range.size)
for i in range(phi_range.size):
    print('t', i, 'phi', phi_range[i])
    t_range[i], error = integrate.quad(t_integral, phi_0, phi_range[i])
for i in range(phi_range.size):
    print('theta', i, 'phi', phi_range[i])
    theta_range[i], error = integrate.quad(theta_integral, phi_0, phi_range[i])

'''print(phi_range)
print(theta_range)
print(t_range)'''

# Computing gradients
theta_grad = np.diff(theta_range) / np.diff(t_range)
phi_grad = np.diff(phi_range) / np.diff(t_range)
theta_grad = np.append(theta_grad, theta_grad[-1])
phi_grad = np.append(phi_grad, phi_grad[-1])

# Visualisation
fig, axs = plt.subplots(2, 2, figsize=(12, 7))

axs[0][0].plot(t_range, theta_range)
axs[0][0].set_yticks(np.arange(np.pi, 3*np.pi, np.pi))
axs[0][0].set_xlabel('Time')
axs[0][0].set_ylabel('Theta')
axs[0][1].plot(t_range, phi_range)
axs[0][1].set_yticks(np.arange(np.pi, 3*np.pi, np.pi))
axs[0][1].set_xlabel('Time')
axs[0][1].set_ylabel('Phi')
axs[1][0].plot(t_range, theta_grad)
# axs[1][0].set_yticks(np.arange(np.pi, 3*np.pi, np.pi))
axs[1][0].set_xlabel('Time')
axs[1][0].set_ylabel('Theta_dot')
axs[1][1].plot(t_range, phi_grad)
# axs[1][1].set_yticks(np.arange(np.pi, 3*np.pi, np.pi))
axs[1][1].set_xlabel('Time')
axs[1][1].set_ylabel('Phi_dot')

plt.show()
plt.close()
