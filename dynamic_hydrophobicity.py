import csv
import os
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.polynomial import polyfit

with open(r'C:\Users\Xin Wenkang\Documents\Scripts\SYPT\data_hydrophobicity\data_hydro.csv', newline='') as f:
    reader = csv.reader(f)
    data = list(reader)

data = data[1:]

omega_range = []
L_range = []
V_range = []
U_range = []
repel = []
crit = []

for l in data:
    omega_range.append(float(l[0])*np.pi/30)
    L_range.append(float(l[1]))
    V_range.append(float(l[2]))
    U_range.append(float(l[-3]))
    repel.append(int(l[-2]))
    crit.append(bool(l[-1]))

omega_range = np.array(omega_range)
L_range = np.array(L_range)
V_range = np.array(V_range)
U_range = np.array(U_range)

crit_index = []
success_index = []
failure_index = []

for i in range(len(crit)):
    if crit[i] == 1:
        crit_index.append(i)
    else:
        pass

for i in range(len(repel)):
    if repel[i] == 1:
        success_index.append(i)
    else:
        failure_index.append(i)

L_success_range = np.array([L_range[i] for i in success_index])
V_success_range = np.array([V_range[i] for i in success_index])
U_success_range = np.array([U_range[i] for i in success_index])

L_failure_range = np.array([L_range[i] for i in failure_index])
V_failure_range = np.array([V_range[i] for i in failure_index])
U_failure_range = np.array([U_range[i] for i in failure_index])

L_crit_range = np.array([L_range[i] for i in crit_index])
V_crit_range = np.array([V_range[i] for i in crit_index])
U_crit_range = np.array([U_range[i] for i in crit_index])

b, m = polyfit(V_crit_range*(L_crit_range**(1/3)), U_crit_range, 1)
x = np.linspace(0, 25, 100)
y = b + m*x

fig, axs = plt.subplots(1, 1, figsize=(12, 7))

axs.plot(V_success_range*(L_success_range**(1/3)), U_success_range, '.', color='limegreen')
axs.plot(V_failure_range*(L_failure_range**(1/3)), U_failure_range, '.', color='brown')
axs.plot(V_crit_range*(L_crit_range**(1/3)), U_crit_range, 'o', color='black')
axs.set_xlabel('$VL^\\frac{1}{3}/m^\\frac{4}{3} \cdot s^{-1}$')
axs.set_ylabel('$U/m \cdot s^{-1}$')
plt.plot(x, y, '-', color='black')
success = plt.fill_between(x,y,color='white')
maxy = plt.ylim()[1]
failure = plt.fill_between(x,y,maxy,color='lightsteelblue')
plt.xlim([0, 25])
plt.ylim([0, y.max()])
plt.legend(handles=[success, failure], labels=['Success', 'Failure'], loc=0)

plt.show()
plt.close()

fig, axs = plt.subplots(1, 1, figsize=(12, 7))

axs.plot(V_success_range, U_success_range, '.', color='limegreen')
axs.plot(V_failure_range, U_failure_range, '.', color='brown')
axs.plot(V_crit_range, U_crit_range, 'o', color='black')
axs.plot(V_crit_range, U_crit_range, color='black')
axs.errorbar(V_crit_range, U_crit_range, xerr=1, color='black')
success = plt.fill_between(V_crit_range,U_crit_range,color='white')
maxy = plt.ylim()[1]
failure = plt.fill_between(V_crit_range,U_crit_range,maxy,color='lightsteelblue')
plt.xlim([V_crit_range.min()-0.5, V_crit_range.max()+0.5])
plt.ylim([0, U_failure_range.max()+0.01])
plt.legend(handles=[success, failure], labels=['Success', 'Failure'], loc=0)
axs.set_xlabel('$V/m^ \cdot s^{-1}$')
axs.set_ylabel('$U/m \cdot s^{-1}$')

plt.savefig('phase_diagram_.png', dpi=300)
plt.show()
plt.close()