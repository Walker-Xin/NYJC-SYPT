import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

# Setting parameters
m = 0.2
m_0 = 0.5
r = 0.1
phi_0 = 30/180 * np.pi
theta_dot_0 = np.pi/1
g = 9.81

I = 0.5 * m_0 * (r**2)
M = (m * (np.sin(phi_0)**2) * (r**2) + I) * theta_dot_0
E = 0.5 * m * (np.sin(phi_0)**2) * (r**2) * (theta_dot_0**2) + 0.5 * I * (theta_dot_0**2) + m * g * r * (1 + np.cos(phi_0))

# Defining integrals
def t_integral(phi):
    c_1 = 2 * E / (m*(r**2))
    c_2 = M**2
    c_3 = (m**2) * (r**4)
    c_4 = I * m * (r**2)
    c_5 = 2 * g / r
    integral = 1 / np.sqrt(c_1 - c_2 / (c_3 * (np.sin(phi)**2) + c_4) - c_5 * (1 + np.cos(phi)))
    return integral

def theta_integral(phi):
    c_1 = 2 * E / (m*r*r)
    c_2 = M**2
    c_3 = (m**2) * (r**4)
    c_4 = I * m * (r**2)
    c_5 = 2 * g / r
    integral_1 = M / ((m * (r**2) * (np.sin(phi)**2)) + I)
    integral_2 = 1 / np.sqrt(c_1 - c_2 / (c_3 * (np.sin(phi)**2) + c_4) - c_5 * (1 + np.cos(phi)))
    integral = integral_1 * integral_2
    return integral

# Setting upper bound
phi_f = 360/180 * np.pi

# Computing integrals
phi_range = np.linspace(phi_0, phi_f, num=100)
t_range = np.zeros(phi_range.size)
theta_range = np.zeros(phi_range.size)
for i in range(phi_range.size): 
    t_range[i], error = integrate.quad(t_integral, phi_0, phi_range[i])
for i in range(phi_range.size): 
    theta_range[i], error = integrate.quad(theta_integral, phi_0, phi_range[i])

print(phi_range, theta_range, t_range)

# Visualisation
fig, axs = plt.subplots(1, 2, figsize=(12, 7))

axs[0].plot(t_range, theta_range)
axs[0].set_yticks(np.arange(np.pi, 3*np.pi, np.pi))
axs[0].set_xlabel('Time')
axs[0].set_ylabel('Theta')
axs[1].plot(t_range, phi_range)
axs[1].set_yticks(np.arange(np.pi, 3*np.pi, np.pi))
axs[1].set_xlabel('Time')
axs[1].set_ylabel('Phi')

plt.show()
plt.close()