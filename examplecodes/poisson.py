from dolfin import *
from dolfin_adjoint import *
import numpy as np
import fenics_optimize as op

m = 0.3    # Target rate of the material amount
p = 5      # Penalty parameter
eps = 1e-3 # Material lower bound
R = 0.1    # Helmholtz filter radius
n = 256    # Resolution

def k(a):
    return eps + (1 - eps) * a ** p

mesh = UnitSquareMesh(n, n)
X = FunctionSpace(mesh, 'CG', 1)
f = interpolate(Constant(1e-2), X)

class Left(SubDomain):
    def inside(self, x, on_boundary):
        gamma = 1/n + eps
        return x[0] == 0.0 and 0.5 - gamma < x[1] < 0.5 + gamma and on_boundary

U = FunctionSpace(mesh, 'CG', 1)
T = Function(U, name='Temperture')
t = TrialFunction(U)
dt = TestFunction(U)
bc = DirichletBC(U, Constant(0.0), Left())

@op.with_derivative([X])
def forward(x):
    rho = op.helmholtzFilter(x[0], X, R=R)
    rho.rename('label', 'control')
    a = inner(grad(t), k(rho)*grad(dt))*dx
    L = f*dt*dx
    A, b = assemble_system(a, L, bc)
    T_s = op.AMGsolver(A, b).solve(T, U, False)
    J = assemble(inner(grad(T), k(rho)*grad(T_s))*dx)
    return J

@op.with_derivative([X])
def constraint(xs):
    rho_bulk = project(Constant(1.0), X)
    rho_0 = assemble(rho_bulk*dx)
    rho_f = assemble(xs[0]*dx)
    rel = rho_f/rho_0
    return rel - m

problemSize = Function(X).vector().size()
x0 = np.ones(problemSize) * m

op.MMAoptimize(problemSize, x0, forward, [constraint], [0], maxeval=100, bounds=[0, 1], rel=1e-20)