from sympy import *
from sympy.physics.mechanics import dynamicsymbols, init_vprinting

t = symbols('t')

psi, theta, phi, x1, x2, x3 = dynamicsymbols('psi, theta, phi, x1, x2, x3')

D_psi, D_theta, D_phi, D_x1, D_x2, D_x3 = dynamicsymbols('D_psi, D_theta, D_phi, D_x1, D_x2, D_x3')

m, g, r, lambdat, lambdaa = symbols('m, g, r, \lambda_t, \lambda_a')

lambdat = m*(r**2)
lambdaa = 0.5*m*(r**2)
I = Matrix([[lambdaa],
            [lambdaa], 
            [lambdat]]) # Intertia tensor

init_vprinting()

# Define 3-1-3 Euler angles
R1 = Matrix([[cos(psi), sin(psi), 0],
            [-sin(psi), cos(psi), 0],
            [0, 0, 1]])

R2 = Matrix([[1, 0, 0],
            [0, cos(theta), sin(theta)],
            [0, -sin(theta), cos(theta)]])

R3 = Matrix([[cos(phi), sin(phi), 0],   
            [-sin(phi), cos(phi), 0],
            [0, 0, 1]])

# Obtain angular velocity
omega = simplify((R3*R2*R1)*Matrix([[0], [0], [psi.diff(t)]]) + (R3*R2)*Matrix([[theta.diff(t)], [0], [0]]) + R3*Matrix([[0], [0], [phi.diff(t)]]))

# Obtain velocity of centre of mass in the spacial coordinate system
R_v = Matrix([[0], [-r], [0]])
v_p = omega.cross(R3*R_v)
v_cm = simplify((R3*R2*R1).T * v_p)

# Solve for the kinematic constraints
Dx1 = -Matrix([[1, 0, 0]]) * v_cm # Instantaneous 'x velocity'
Dx2 = -Matrix([[0, 1, 0]]) * v_cm # Instantaneous 'y velocity'
# Dx3 = n_z.dot(v_cm)

T = 0.5*m*simplify(Dx1**2 + Dx2**2)[0] + 0.5*I.dot(diag(omega[0], omega[1], omega[2])*omega) # Kinematic energy
V = m*g*x3

L = T - V

pprint(L)