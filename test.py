import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import scipy.optimize as optimize

omega = np.pi/0.25 # Angular velocity of the loop in radian per second
g = 9.81  # Gravitational accelertaion in meter per second^2
r = 0.1  # Radius of the loop in meter
mu = 0.2

A = omega**2
B = g / r

def f_plus(phi):
    return B - A*np.cos(phi) + mu*(B*(1/np.sin(phi)-1/np.tan(phi)) - A*np.sin(phi))

def f_minus(phi):
    return B - A*np.cos(phi) + mu*(B*(1/np.sin(phi)-1/np.tan(phi)) - A*np.sin(phi))

phi_sample = np.linspace(0.01, 10, 100)

f = np.zeros(np.size(phi_sample))

for i in range(len(phi_sample)):
    f[i] = f_minus(phi_sample[i])

plt.plot(phi_sample, f)
plt.plot(phi_sample, np.zeros(np.size(phi_sample)))
plt.show()

phi_e = optimize.root_scalar(f_plus, bracket=[0.01, 1.5], method='brentq')
print(phi_e.root)