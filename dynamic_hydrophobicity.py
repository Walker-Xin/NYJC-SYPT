import csv
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.polynomial import polyfit

r_s = 1.0
r_l = 1.6

with open(r'C:\Users\Xin Wenkang\OneDrive\Files\Documents\SYPT\Dynamic Hydrophobicity\data\data_hydro_large.csv', newline='') as f:
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

L_success_range_l = np.array([L_range[i] for i in success_index])
V_success_range_l = np.array([V_range[i] for i in success_index])
U_success_range_l = np.array([U_range[i] for i in success_index])

L_failure_range_l = np.array([L_range[i] for i in failure_index])
V_failure_range_l = np.array([V_range[i] for i in failure_index])
U_failure_range_l = np.array([U_range[i] for i in failure_index])

L_crit_range_l = np.array([L_range[i] for i in crit_index])
V_crit_range_l = np.array([V_range[i] for i in crit_index])
U_crit_range_l = np.array([U_range[i] for i in crit_index])

b, m = polyfit(V_crit_range_l*(L_crit_range_l**(1/3)), U_crit_range_l**(4/3), 1)
m_l = 0.02159415241010928
x_l = np.linspace(0, 25, 100)
y_l = m_l*x_l
print(m_l)

fig, axs = plt.subplots(1, 1, figsize=(8, 6))

axs.plot(V_success_range_l*(L_success_range_l**(1/3)), U_success_range_l**(4/3), '.', color='limegreen')
axs.plot(V_failure_range_l*(L_failure_range_l**(1/3)), U_failure_range_l**(4/3), '.', color='brown')
axs.plot(V_crit_range_l*(L_crit_range_l**(1/3)), U_crit_range_l**(4/3), 'o', color='black')
axs.errorbar(V_crit_range_l*(L_crit_range_l**(1/3)), U_crit_range_l**(4/3), yerr=0.01, fmt='o', capsize=3, color='black')
axs.set_xlabel('$VL^\\frac{1}{3}$')
axs.set_ylabel('$U^\\frac{4}{3}$')
plt.plot(x_l, y_l, '-', color='black')
success = plt.fill_between(x_l,y_l,color='white')
maxy = plt.ylim()[1]
failure = plt.fill_between(x_l,y_l,maxy,color='lightsteelblue')
plt.xlim([0, 25])
plt.ylim([0, y_l.max()])
plt.legend(handles=[success, failure], labels=['Success', 'Failure'], loc=0)
plt.text(0.85, 0.01, r'$r$ = {} mm'.format(r_l), transform=plt.gca().transAxes)

plt.show()
plt.close()

fig, axs = plt.subplots(1, 1, figsize=(8, 6))

axs.plot(V_success_range_l, U_success_range_l, '.', color='limegreen')
axs.plot(V_failure_range_l, U_failure_range_l, '.', color='brown')
axs.plot(V_crit_range_l, U_crit_range_l, 'o', color='black')
axs.plot(V_crit_range_l, U_crit_range_l, color='black')
axs.errorbar(V_crit_range_l, U_crit_range_l, yerr=0.01, fmt='o', capsize=3, color='black')
success = plt.fill_between(V_crit_range_l,U_crit_range_l,color='white')
maxy = plt.ylim()[1]
failure = plt.fill_between(V_crit_range_l,U_crit_range_l,maxy,color='lightsteelblue')
plt.xlim([V_crit_range_l.min()-0.5, V_crit_range_l.max()+0.5])
plt.ylim([0, U_failure_range_l.max()+0.01])
plt.legend(handles=[success, failure], labels=['Success', 'Failure'], loc=0)
axs.set_xlabel('$V/m^ \cdot s^{-1}$')
axs.set_ylabel('$U/m \cdot s^{-1}$')
plt.text(0.85, 0.01, r'$r$ = {} mm'.format(r_l), transform=plt.gca().transAxes)

# plt.savefig('phase_diagram_.png', dpi=300)
plt.show()
plt.close()



with open(r'C:\Users\Xin Wenkang\OneDrive\Files\Documents\SYPT\Dynamic Hydrophobicity\data\data_hydro_small.csv', newline='') as f:
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

L_success_range_s = np.array([L_range[i] for i in success_index])
V_success_range_s = np.array([V_range[i] for i in success_index])
U_success_range_s = np.array([U_range[i] for i in success_index])

L_failure_range_s = np.array([L_range[i] for i in failure_index])
V_failure_range_s = np.array([V_range[i] for i in failure_index])
U_failure_range_s = np.array([U_range[i] for i in failure_index])

L_crit_range_s = np.array([L_range[i] for i in crit_index])
V_crit_range_s = np.array([V_range[i] for i in crit_index])
U_crit_range_s = np.array([U_range[i] for i in crit_index])

b, m = polyfit(V_crit_range_s*(L_crit_range_s**(1/3)), U_crit_range_s**(4/3), 1)
m_s = 1.1*0.0257613160788931
x_s = np.linspace(0, 25, 100)
y_s = m_s*x_s
print(m_s)

fig, axs = plt.subplots(1, 1, figsize=(8, 6))

axs.plot(V_success_range_s*(L_success_range_s**(1/3)), U_success_range_s**(4/3), '.', color='limegreen')
axs.plot(V_failure_range_s*(L_failure_range_s**(1/3)), U_failure_range_s**(4/3), '.', color='brown')
axs.plot(V_crit_range_s*(L_crit_range_s**(1/3)), U_crit_range_s**(4/3), 'o', color='black')
axs.errorbar(V_crit_range_s*(L_crit_range_s**(1/3)), U_crit_range_s**(4/3), yerr=0.01, fmt='o', capsize=3, color='black')
axs.set_xlabel('$VL^\\frac{1}{3}$')
axs.set_ylabel('$U^\\frac{4}{3}$')
plt.plot(x_s, y_s, '-', color='black')
success = plt.fill_between(x_s,y_s,color='white')
maxy = plt.ylim()[1]
failure = plt.fill_between(x_s,y_s,maxy,color='lightsteelblue')
plt.xlim([0, 25])
plt.ylim([0, y_s.max()])
plt.legend(handles=[success, failure], labels=['Success', 'Failure'], loc=0)
plt.text(0.85, 0.01, r'$r$ = {} mm'.format(r_s), transform=plt.gca().transAxes)

plt.show()
plt.close()

fig, axs = plt.subplots(1, 1, figsize=(8, 6))

axs.plot(V_success_range_s, U_success_range_s, '.', color='limegreen')
axs.plot(V_failure_range_s, U_failure_range_s, '.', color='brown')
axs.plot(V_crit_range_s, U_crit_range_s, 'o', color='black')
axs.plot(V_crit_range_s, U_crit_range_s, color='black')
axs.errorbar(V_crit_range_s, U_crit_range_s, yerr=0.01, fmt='o', capsize=3, color='black')
success = plt.fill_between(V_crit_range_s,U_crit_range_s,color='white')
maxy = plt.ylim()[1]
failure = plt.fill_between(V_crit_range_s,U_crit_range_s,maxy,color='lightsteelblue')
plt.xlim([V_crit_range_s.min()-0.5, V_crit_range_s.max()+0.5])
plt.ylim([0, U_failure_range_s.max()+0.01])
plt.legend(handles=[success, failure], labels=['Success', 'Failure'], loc=0)
axs.set_xlabel('$V/m^ \cdot s^{-1}$')
axs.set_ylabel('$U/m \cdot s^{-1}$')
plt.text(0.85, 0.01, r'$r$ = {} mm'.format(r_s), transform=plt.gca().transAxes)

# plt.savefig('phase_diagram_.png', dpi=300)
plt.show()
plt.close()



fig, axs = plt.subplots(1, 1, figsize=(8, 6))

axs.plot(V_crit_range_l*(L_crit_range_l**(1/3)), U_crit_range_l**(4/3), 'o', color='tab:blue', label='large')

axs.plot(V_crit_range_s*(L_crit_range_s**(1/3)), U_crit_range_s**(4/3), 'o', color='tab:red', label='small')

plt.plot(x_l, y_l, '-', color='tab:blue')
plt.plot(x_s, y_s, '-', color='tab:red')

axs.set_xlabel('$VL^\\frac{1}{3}$')
axs.set_ylabel('$U^\\frac{4}{3}$')

plt.legend(labels=['$r$ = 1.6mm', '$r$ = 1.0mm'])

plt.xlim([0, 27])
plt.ylim([0, y_s.max()])

plt.show()
plt.close()
print(m_s/m_l)



with open(r'C:\Users\Xin Wenkang\OneDrive\Files\Documents\SYPT\Dynamic Hydrophobicity\data\motion_large.csv', newline='') as f:
    reader = csv.reader(f)
    data = list(reader)

frame_range = []
x_range = []
y_range = []
i = 0

for l in data:
    frame_range.append(i)
    x_range.append(float(l[1]))
    y_range.append(float(l[2]))
    i += 1

frame_range = np.array(frame_range)
t_range_l = frame_range*(1/800)
x_range_l = np.array(x_range)
y_range_l = np.array(y_range)

fig, axs = plt.subplots(1, 1, figsize=(8, 6))

axs.plot(x_range_l, y_range_l, '.')
axs.set_xlabel('x/cm')
axs.set_ylabel('y/cm')
plt.text(0.85, 0.01, r'$r$ = {} mm'.format(r_l), transform=plt.gca().transAxes)

plt.show()
plt.close()

b, m = polyfit(np.extract(x_range_l>0, t_range_l), np.extract(x_range_l>0, x_range_l), 1)
x = np.linspace(np.extract(x_range_l>0, t_range_l).min(), np.extract(x_range_l>0, t_range_l).max(), 100)
y = b + m*x
print(m)

fig, axs = plt.subplots(1, 1, figsize=(8, 6))

axs.plot(t_range_l, x_range_l, '.')
plt.plot(x, y, '-', color='tab:red')
axs.set_xlabel('t/s')
axs.set_ylabel('x/cm')
plt.text(0.80, 0.01, r'$r$ = {} mm'.format(r_l), transform=plt.gca().transAxes)
plt.text(0.80, 0.06, r'$\Delta u_x$ = {} cm/s'.format(round(m, 1)), transform=plt.gca().transAxes)

plt.show()
plt.close()



with open(r'C:\Users\Xin Wenkang\OneDrive\Files\Documents\SYPT\Dynamic Hydrophobicity\data\motion_small.csv', newline='') as f:
    reader = csv.reader(f)
    data = list(reader)

frame_range = []
x_range = []
y_range = []
i = 0

for l in data:
    frame_range.append(i)
    x_range.append(float(l[1]))
    y_range.append(float(l[2]))
    i += 1

frame_range = np.array(frame_range)
t_range_s = frame_range*(1/800)
x_range_s = np.array(x_range)
y_range_s = np.array(y_range)

fig, axs = plt.subplots(1, 1, figsize=(8, 6))

axs.plot(x_range_s, y_range_s, '.')
axs.set_xlabel('x/cm')
axs.set_ylabel('y/cm')
plt.text(0.85, 0.01, r'$r$ = {} mm'.format(r_s), transform=plt.gca().transAxes)

plt.show()
plt.close()

b, m = polyfit(np.extract(x_range_s>0, t_range_s), np.extract(x_range_s>0, x_range_s), 1)
x = np.linspace(np.extract(x_range_s>0, t_range_s).min(), np.extract(x_range_s>0, t_range_s).max(), 100)
y = b + m*x
print(m)

fig, axs = plt.subplots(1, 1, figsize=(8, 6))

axs.plot(t_range_s, x_range_s, '.')
plt.plot(x, y, '-', color='tab:red')
axs.set_xlabel('t/s')
axs.set_ylabel('x/cm')
plt.text(0.80, 0.01, r'$r$ = {} mm'.format(r_s), transform=plt.gca().transAxes)
plt.text(0.80, 0.06, r'$\Delta u_x$ = {} cm/s'.format(round(m, 1)), transform=plt.gca().transAxes)

plt.show()
plt.close()



with open(r'C:\Users\Xin Wenkang\OneDrive\Files\Documents\SYPT\Dynamic Hydrophobicity\data\speed_incre_large.csv', newline='') as f:
    reader = csv.reader(f)
    data = list(reader)

data = data[1:]

V_l = []
ux_l = []

for l in data:
    V_l.append(float(l[2]))
    ux_l.append(float(l[-1]))

V_l = np.array(V_l)
ux_l = np.array(ux_l)

with open(r'C:\Users\Xin Wenkang\OneDrive\Files\Documents\SYPT\Dynamic Hydrophobicity\data\speed_incre_small.csv', newline='') as f:
    reader = csv.reader(f)
    data = list(reader)

data = data[1:]

V_s = []
ux_s = []

for l in data:
    V_s.append(float(l[2]))
    ux_s.append(float(l[-1]))

V_s = np.array(V_s)
ux_s = np.array(ux_s)

a_l = 0.8*0.000643
x_l = np.linspace(0, V_l.max()+5, 100)
y_l = a_l*(x_l**2)

a_s = 1.25*0.0011628
x_s = np.linspace(0, V_s.max()+5, 100)
y_s = a_s*(x_s**2)

fig, axs = plt.subplots(1, 1, figsize=(8, 6))

axs.plot(V_s, ux_s, 'o', color='tab:red')
axs.plot(V_l, ux_l, 'o', color='tab:blue')
axs.errorbar(V_s, ux_s, yerr=0.1, fmt='o', capsize=3, color='tab:red')
axs.errorbar(V_l, ux_l, yerr=0.1, fmt='o', capsize=3, color='tab:blue')
axs.plot(x_s, y_s, color='tab:red')
axs.plot(x_l, y_l, color='tab:blue')
# plt.xlim([0, V_l.max()+3])
# plt.ylim([0, ux_l.max()+0.1])
axs.set_xlabel('$V/m^{2} \cdot s^{-2}$')
axs.set_ylabel('$\Delta u_x/m \cdot s^{-1}$')
plt.legend(labels=['$r$ = 1.0mm', '$r$ = 1.6mm'])
# plt.text(0.80, 0.01, r'$a_l$ = {}'.format(round(a_l,6)), transform=plt.gca().transAxes)
# plt.text(0.80, 0.06, r'$a_s$ = {}'.format(round(a_s,6)), transform=plt.gca().transAxes)

plt.show()
plt.close()